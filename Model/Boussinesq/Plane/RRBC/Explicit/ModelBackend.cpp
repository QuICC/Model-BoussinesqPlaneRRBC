/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Plane/RRBC/Explicit/ModelBackend.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitBoundaryValue.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Heating.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Lapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Lapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Lapl2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Id.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

namespace Explicit {

ModelBackend::ModelBackend() :
    IRRBCBackend(),
#ifdef QUICC_TRANSFORM_CHEBYSHEV_TRUNCATE_QI
    mcTruncateQI(true)
#else
    mcTruncateQI(false)
#endif // QUICC_TRANSFORM_CHEBYSHEV_TRUNCATE_QI
{}

bool ModelBackend::isComplex(const SpectralFieldId& fId) const
{
   return false;
}

ModelBackend::SpectralFieldIds ModelBackend::implicitFields(
   const SpectralFieldId& fId) const
{
   auto velTor = std::make_pair(PhysicalNames::Velocity::id(),
      FieldComponents::Spectral::TOR);
   auto velPol = std::make_pair(PhysicalNames::Velocity::id(),
      FieldComponents::Spectral::POL);
   auto temp = std::make_pair(PhysicalNames::Temperature::id(),
      FieldComponents::Spectral::SCALAR);
   SpectralFieldIds fields = {velTor, velPol, temp};

   return fields;
}

ModelBackend::SpectralFieldIds ModelBackend::explicitNonlinearFields(
   const SpectralFieldId& fId) const
{
   SpectralFieldIds fields;
   if (fId == std::make_pair(PhysicalNames::Temperature::id(),
                 FieldComponents::Spectral::SCALAR))
   {
      fields.push_back(std::make_pair(PhysicalNames::Temperature::id(),
         FieldComponents::Spectral::SCALAR));
   }

   return fields;
}

void ModelBackend::equationInfo(EquationInfo& info, const SpectralFieldId& fId,
   const Resolution& res) const
{
   // Operators are real
   info.isComplex = this->isComplex(fId);

   // Splitting 4th poloidal equation into two systems
   if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                 FieldComponents::Spectral::POL))
   {
      info.isSplitEquation = this->useSplitEquation();
   }
   else
   {
      info.isSplitEquation = false;
   }

   // Implicit coupled fields
   info.im = this->implicitFields(fId);

   // Explicit linear terms
   info.exL.clear();

   // Explicit nonlinear terms
   info.exNL = this->explicitNonlinearFields(fId);

   // Explicit nextstep terms
   info.exNS.clear();

   // Index mode
   info.indexMode = static_cast<int>(Equations::CouplingIndexType::MODE);
}

void ModelBackend::operatorInfo(OperatorInfo& info, const SpectralFieldId& fId,
   const Resolution& res, const Equations::Tools::ICoupling& coupling,
   const BcMap& bcs) const
{
   // Loop overall matrices/eigs
   for (int idx = 0; idx < info.tauN.size(); ++idx)
   {
      auto eigs = coupling.getIndexes(res, idx);

      int tN, gN, rhs;
      ArrayI shift(3);

      this->blockInfo(tN, gN, shift, rhs, fId, res, idx, bcs);

      info.tauN(idx) = tN;
      info.galN(idx) = gN;
      info.galShift.row(idx) = shift;
      info.rhsCols(idx) = rhs;

      // Compute system size
      int sN = 0;
      for (auto f: this->implicitFields(fId))
      {
         this->blockInfo(tN, gN, shift, rhs, f, res, idx, bcs);
         sN += gN;
      }

      if (sN == 0)
      {
         sN = info.galN(idx);
      }

      info.sysN(idx) = sN;
   }
}

std::vector<details::BlockDescription> ModelBackend::implicitBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplitOperator) const
{
   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->zi = nds.find(NonDimensional::Lower1d::id())->second->value();
      opts->zo = nds.find(NonDimensional::Upper1d::id())->second->value();
      opts->k1 = eigs.at(0);
      opts->k2 = eigs.at(1);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = isSplitOperator;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR))
   {
      if (rowId == colId)
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int k1,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            if (o.k1 == 0 && o.k2 == 0)
            {
               SparseSM::Chebyshev::LinearMap::Id spasm(nNr, nNc, o.zi, o.zo, 2,
                  0);
               bMat = spasm.mat();
            }
            else
            {
               auto laplh = -(o.k1 * o.k1 + o.k2 * o.k2);
               SparseSM::Chebyshev::LinearMap::I2Lapl spasm(nNr, nNc, o.zi,
                  o.zo, o.k1, o.k2);
               bMat = laplh * spasm.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else if (colId == std::make_pair(PhysicalNames::Velocity::id(),
                           FieldComponents::Spectral::POL))
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int k1,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            auto E = nds.find(NonDimensional::Ekman::id())->second->value();

            if (o.k1 == 0 && o.k2 == 0)
            {
               auto c = 1. / E;
               SparseSM::Chebyshev::LinearMap::I2 spasm(nNr, nNc, o.zi, o.zo);
               bMat = c * spasm.mat();
            }
            else
            {
               auto c = -(o.k1 * o.k1 + o.k2 * o.k2) / E;
               SparseSM::Chebyshev::LinearMap::I2D1 i2d1(nNr, nNc, o.zi, o.zo);
               bMat = c * i2d1.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                        FieldComponents::Spectral::POL))
   {
      if (rowId == colId)
      {
         // Real part of block
         auto realOp = [](const int nNr, const int nNc, const int k1,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);
            auto laplh = -(o.k1 * o.k1 + o.k2 * o.k2);
            if (o.useSplitEquation)
            {
               if (o.isSplitOperator)
               {
                  SparseSM::Chebyshev::LinearMap::I2Lapl spasm(nNr, nNc, o.zi,
                     o.zo, o.k1, o.k2);
                  bMat = laplh * spasm.mat();
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I2Lapl spasm(nNr, nNc, o.zi,
                     o.zo, o.k1, o.k2);
                  bMat = laplh * spasm.mat();
               }
            }
            else
            {
               if (o.k1 == 0 && o.k2 == 0)
               {
                  SparseSM::Chebyshev::LinearMap::Id spasm(nNr, nNc, o.zi, o.zo,
                     2, 0);
                  bMat = spasm.mat();
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I4Lapl2 spasm(nNr, nNc, o.zi,
                     o.zo, o.k1, o.k2);
                  bMat = laplh * spasm.mat();
               }
            }

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else if (colId == std::make_pair(PhysicalNames::Velocity::id(),
                           FieldComponents::Spectral::TOR))
      {
         // Real part of block
         auto realOp = [](const int nNr, const int nNc, const int k1,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);
            auto E = nds.find(NonDimensional::Ekman::id())->second->value();

            if (o.k1 == 0 && o.k2 == 0)
            {
               auto c = -1. / E;
               SparseSM::Chebyshev::LinearMap::I2 i2(nNr, nNc, o.zi, o.zo);
               bMat = c * i2.mat();
            }
            else
            {
               auto c = (o.k1 * o.k1 + o.k2 * o.k2) / E;
               SparseSM::Chebyshev::LinearMap::I4D1 i4d1(nNr, nNc, o.zi, o.zo);
               bMat = c * i4d1.mat();
            }

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else if (colId == std::make_pair(PhysicalNames::Temperature::id(),
                           FieldComponents::Spectral::SCALAR))
      {
         // Real part of block
         auto realOp = [](const int nNr, const int nNc, const int k1,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);
            auto Ra = nds.find(NonDimensional::Rayleigh::id())->second->value();
            auto Pr = nds.find(NonDimensional::Prandtl::id())->second->value();

            auto laplh = -(o.k1 * o.k1 + o.k2 * o.k2);
            SparseSM::Chebyshev::LinearMap::I4 spasm(nNr, nNc, o.zi, o.zo);
            bMat = -(Ra / Pr) * laplh * spasm.mat();

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
                        FieldComponents::Spectral::SCALAR))
   {
      if (rowId == colId)
      {
         // Creat real part of block
         auto realOp = [](const int nNr, const int nNc, const int k1,
                          std::shared_ptr<details::BlockOptions> opts,
                          const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto Pr =
               nds.find(NonDimensional::Prandtl::id())->second->value();

            if (o.k1 == 0 && o.k2 == 0)
            {
               SparseSM::Chebyshev::LinearMap::Id spasm(nNr, nNc, o.zi, o.zo, 2,
                  0);
               bMat = (1.0 / Pr) * spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2Lapl spasm(nNr, nNc, o.zi,
                  o.zo, o.k1, o.k2);
               bMat = (1.0 / Pr) * spasm.mat();
            }

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
   }
   else
   {
      throw std::logic_error("Equations are not setup properly [(" +
                             PhysicalNames::Coordinator::tag(rowId.first) +
                             ", " + std::to_string(rowId.second) + "),(" +
                             PhysicalNames::Coordinator::tag(colId.first) +
                             ", " + std::to_string(colId.second) + ")]");
   }

   return descr;
}

std::vector<details::BlockDescription> ModelBackend::timeBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(rowId == colId);
   auto fieldId = rowId;

   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->zi = nds.find(NonDimensional::Lower1d::id())->second->value();
      opts->zo = nds.find(NonDimensional::Upper1d::id())->second->value();
      opts->k1 = eigs.at(0);
      opts->k2 = eigs.at(1);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = false;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                     FieldComponents::Spectral::TOR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int k1,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (o.k1 == 0 && o.k2 == 0)
         {
            SparseSM::Chebyshev::LinearMap::I2 spasm(nNr, nNc, o.zi, o.zo);
            bMat = spasm.mat();
         }
         else
         {
            auto laplh = -(o.k1 * o.k1 + o.k2 * o.k2);
            SparseSM::Chebyshev::LinearMap::I2 spasm(nNr, nNc, o.zi, o.zo);
            bMat = laplh * spasm.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                          FieldComponents::Spectral::POL))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int k1,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat;
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         if (o.useSplitEquation)
         {
            SparseSM::Chebyshev::LinearMap::I2 spasm(nNr, nNc, o.zi, o.zo);
            bMat = spasm.mat();
         }
         else
         {
            if (o.k1 == 0 && o.k2 == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2 spasm(nNr, nNc, o.zi, o.zo);
               bMat = spasm.mat();
            }
            else
            {
               auto laplh = -(o.k1 * o.k1 + o.k2 * o.k2);
               SparseSM::Chebyshev::LinearMap::I4Lapl spasm(nNr, nNc, o.zi,
                  o.zo, o.k1, o.k2);

#if 0
               // Correct Laplacian for 4th order system according to:
               // McFadden,Murray,Boisvert,
               // Elimination of Spurious Eigenvalues in the
               // Chebyshev Tau Spectral Method,
               // JCP 91, 228-239 (1990)
               // We simply drop the last two column
               if (o.bcId == Bc::Name::NoSlip::id())
               {
                  SparseSM::Chebyshev::LinearMap::Id qid(nNr, nNc, o.zi, o.zo, -2);
                  bMat = laplh * spasm.mat() * qid.mat();
               }
               else
               {
                  bMat = laplh * spasm.mat();
               }
#else
               bMat = laplh * spasm.mat();
#endif
            }
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else if (fieldId == std::make_pair(PhysicalNames::Temperature::id(),
                          FieldComponents::Spectral::SCALAR))
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int k1,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         SparseMatrix bMat;
         if (o.k1 == 0 && o.k2 == 0)
         {
            SparseSM::Chebyshev::LinearMap::I2 spasm(nNr, nNc, o.zi, o.zo);
            bMat = spasm.mat();
         }
         else
         {
            SparseSM::Chebyshev::LinearMap::I2 spasm(nNr, nNc, o.zi, o.zo);
            bMat = spasm.mat();
         }

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }

   return descr;
}

std::vector<details::BlockDescription> ModelBackend::boundaryBlockBuilder(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplit) const
{
   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->zi = nds.find(NonDimensional::Lower1d::id())->second->value();
      opts->zo = nds.find(NonDimensional::Upper1d::id())->second->value();
      opts->k1 = eigs.at(0);
      opts->k2 = eigs.at(1);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = isSplit;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (rowId == colId)
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int k1,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         SparseMatrix bMat(nNr, nNc);

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }

   return descr;
}

std::vector<details::BlockDescription>
ModelBackend::splitBoundaryValueBlockBuilder(const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(rowId == colId);
   auto fieldId = rowId;

   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->zi = nds.find(NonDimensional::Lower1d::id())->second->value();
      opts->zo = nds.find(NonDimensional::Upper1d::id())->second->value();
      opts->k1 = eigs.at(0);
      opts->k2 = eigs.at(1);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = false;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                     FieldComponents::Spectral::POL))
   {
      // Boundary value operator
      auto bcValOp = [](const int nNr, const int nNc, const int k1,
                        std::shared_ptr<details::BlockOptions> opts,
                        const NonDimensional::NdMap& nds)
      {
         assert(nNr == nNc);

         SparseMatrix bMat(nNr, 2);

         Eigen::Triplet<MHDFloat> valTop = {0, 0, 1.0};
         Eigen::Triplet<MHDFloat> valBot = {1, 1, 1.0};
         std::vector<Eigen::Triplet<MHDFloat>> triplets = {valTop, valBot};
         bMat.setFromTriplets(triplets.begin(), triplets.end());

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = bcValOp;
      d.imagOp = bcValOp;
   }

   return descr;
}

void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix,
   const std::size_t opId,
   const Equations::CouplingInformation::FieldId_range imRange,
   const int matIdx, const std::size_t bcType, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(eigs.size() == 2);
   int k = res.cpu()->dim(Dimensions::Transform::SPECTRAL)->mode(matIdx)(0);

   // Time operator
   if (opId == ModelOperator::Time::id())
   {
      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         auto colId = rowId;
         const auto& fields = this->implicitFields(rowId);
         auto descr = timeBlockBuilder(rowId, colId, res, eigs, bcs, nds);
         buildBlock(rModelMatrix, descr, rowId, colId, fields, matIdx, bcType,
            res, k, k, bcs, nds, false);
      }
   }
   // Linear operator
   else if (opId == ModelOperator::ImplicitLinear::id() ||
            opId == ModelOperator::SplitImplicitLinear::id())
   {
      bool isSplit = (opId == ModelOperator::SplitImplicitLinear::id());

      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr =
               implicitBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
            buildBlock(rModelMatrix, descr, rowId, colId, fields, matIdx,
               bcType, res, k, k, bcs, nds, isSplit);
         }
      }
   }
   // Boundary operator
   else if (opId == ModelOperator::Boundary::id() ||
            opId == ModelOperator::SplitBoundary::id())
   {
      bool isSplit = (opId == ModelOperator::SplitBoundary::id());

      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr =
               boundaryBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
            buildBlock(rModelMatrix, descr, rowId, colId, fields, matIdx,
               bcType, res, k, k, bcs, nds, isSplit);
         }
      }
   }
   // Split equation boundary value
   else if (opId == ModelOperator::SplitBoundaryValue::id())
   {
      for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
      {
         auto rowId = *pRowId;
         const auto& fields = this->implicitFields(rowId);
         for (auto pColId = imRange.first; pColId != imRange.second; pColId++)
         {
            auto colId = *pColId;
            auto descr = splitBoundaryValueBlockBuilder(rowId, colId, res, eigs,
               bcs, nds);
            buildFixedBlock(rModelMatrix, 1, true, descr, rowId, colId, fields,
               matIdx, bcType, res, k, k, bcs, nds, false);
         }
      }
   }
   else
   {
      throw std::logic_error("Requested operator type is not implemented");
   }
}

void ModelBackend::galerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& fieldId, const int matIdx, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(eigs.size() == 2);
   int k = res.cpu()->dim(Dimensions::Transform::SPECTRAL)->mode(matIdx)(0);

   this->stencil(mat, fieldId, k, res, makeSquare, bcs, nds);
}

void ModelBackend::explicitBlock(DecoupledZSparse& decMat,
   const SpectralFieldId& rowId, const std::size_t opId,
   const SpectralFieldId colId, const int matIdx, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(eigs.size() == 2);
   int k = res.cpu()->dim(Dimensions::Transform::SPECTRAL)->mode(matIdx)(0);

   auto bcType = ModelOperatorBoundary::SolverNoTau::id();

   // Explicit linear operator
   if (opId == ModelOperator::ExplicitLinear::id())
   {
      throw std::logic_error("There are no explicit linear operators");
   }
   // Explicit nonlinear operator
   else if (opId == ModelOperator::ExplicitNonlinear::id())
   {
      const auto& fields = this->explicitNonlinearFields(rowId);
      auto descr =
         explicitNonlinearBlockBuilder(rowId, colId, res, eigs, bcs, nds);
      buildBlock(decMat, descr, rowId, colId, fields, matIdx, bcType, res, k, k,
         bcs, nds, false, true);
   }
   // Explicit nextstep operator
   else if (opId == ModelOperator::ExplicitNextstep::id())
   {
      throw std::logic_error("There are no explicit nextstep operators");
   }
}

std::vector<details::BlockDescription>
ModelBackend::explicitNonlinearBlockBuilder(const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   assert(eigs.size() == 2);

   std::vector<details::BlockDescription> descr;

   // Create description with common options
   auto getDescription = [&]() -> details::BlockDescription&
   {
      descr.push_back({});
      auto& d = descr.back();
      auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
      opts->zi = nds.find(NonDimensional::Lower1d::id())->second->value();
      opts->zo = nds.find(NonDimensional::Upper1d::id())->second->value();
      opts->k1 = eigs.at(0);
      opts->k2 = eigs.at(1);
      opts->bcId = bcs.find(colId.first)->second;
      opts->truncateQI = this->mcTruncateQI;
      opts->isSplitOperator = false;
      opts->useSplitEquation = this->useSplitEquation();
      d.opts = opts;

      return d;
   };

   if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
                   FieldComponents::Spectral::SCALAR) &&
       rowId == colId)
   {
      // Real part of operator
      auto realOp = [](const int nNr, const int nNc, const int k1,
                       std::shared_ptr<details::BlockOptions> opts,
                       const NonDimensional::NdMap& nds)
      {
         SparseMatrix bMat(nNr, nNc);

         auto& o =
            *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

         SparseSM::Chebyshev::LinearMap::I2 spasm(nNr, nNc, o.zi, o.zo);
         bMat = spasm.mat();

         return bMat;
      };

      // Create block diagonal operator
      auto& d = getDescription();
      d.nRowShift = 0;
      d.nColShift = 0;
      d.realOp = realOp;
      d.imagOp = nullptr;
   }
   else
   {
      throw std::logic_error("There are no explicit nonlinear operators");
   }

   return descr;
}


} // namespace Explicit
} // namespace RRBC
} // namespace Plane
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
