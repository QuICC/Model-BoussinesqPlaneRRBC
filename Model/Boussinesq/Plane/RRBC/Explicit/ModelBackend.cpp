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
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitBoundaryValue.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Heating.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Id.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Lapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Lapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Lapl2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4D1.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

namespace Explicit {

   ModelBackend::ModelBackend()
      : IRRBCBackend(),
#ifdef QUICC_TRANSFORM_CHEBYSHEV_TRUNCATE_QI
      mcTruncateQI(true)
#else
      mcTruncateQI(false)
#endif // QUICC_TRANSFORM_CHEBYSHEV_TRUNCATE_QI
   {
   }

   ModelBackend::SpectralFieldIds ModelBackend::implicitFields(const SpectralFieldId& fId) const
   {
      auto velTor = std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR);
      auto velPol = std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR);
      auto temp = std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR);
      SpectralFieldIds fields = {velTor, velPol, temp};

      return fields;
   }

   ModelBackend::SpectralFieldIds ModelBackend::explicitNonlinearFields(const SpectralFieldId& fId) const
   {
      SpectralFieldIds fields;
      if(fId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
      {
         fields.push_back(std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR));
      }

      return fields;
   }

   void ModelBackend::equationInfo(EquationInfo& info, const SpectralFieldId& fId, const Resolution& res) const
   {
      // Operators are real
      info.isComplex = false;

      // Splitting 4th poloidal equation into two systems
      if(fId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL))
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
      info.indexMode = static_cast<int>(Equations::CouplingIndexType::SLOWEST_MULTI_RHS);
   }

   void ModelBackend::blockSize(int& tN, int& gN, ArrayI& shift, int& rhs, const SpectralFieldId& fId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs) const
   {
      assert(eigs.size() == 2);
      int k1 = eigs.at(0);
      int k2 = eigs.at(1);

      this->blockInfo(tN, gN, shift, rhs, fId, res, k1, k2, bcs);
   }

   void ModelBackend::operatorInfo(OperatorInfo& info, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const
   {
      // Loop overall matrices/eigs
      for(int idx = 0; idx < info.tauN.size(); ++idx)
      {
         auto eigs = coupling.getIndexes(res, idx);

         int tN, gN, rhs;
         ArrayI shift(3);

         this->blockSize(tN, gN, shift, rhs, fId, res, eigs, bcs);

         info.tauN(idx) = tN;
         info.galN(idx) = gN;
         info.galShift.row(idx) = shift;
         info.rhsCols(idx) = rhs;

         // Compute system size
         int sN = 0;
         for(auto f: this->implicitFields(fId))
         {
            this->blockSize(tN, gN, shift, rhs, f, res, eigs, bcs);
            sN += gN;
         }

         if(sN == 0)
         {
            sN = info.galN(idx);
         }

         info.sysN(idx) = sN;
      }
   }

   void ModelBackend::implicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds, const bool isSplitOperator) const
   {
      assert(eigs.size() == 2);
      int k1 = eigs.at(0);
      int k2 = eigs.at(1);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, k1)(0);

      auto zi = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto zo = nds.find(NonDimensional::Upper1d::id())->second->value();

      if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
      {
         if(rowId == colId)
         {
            if(k1 == 0 && k2 == 0)
            {
               SparseSM::Chebyshev::LinearMap::Id spasm(nN, nN, zi, zo, -2, 0);
               decMat.real() = spasm.mat();
            }
            else
            {
               auto laplh = -(k1*k1 + k2*k2);
               SparseSM::Chebyshev::LinearMap::I2Lapl i2lapl(nN, nN, zi, zo, k1, k2);
               decMat.real() = laplh*i2lapl.mat();
            }
         }
         else if(colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
         {
            if(k1 == 0 && k2 == 0)
            {
               auto E = nds.find(NonDimensional::Ekman::id())->second->value();

               auto c = 1./E;
               SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, zi, zo);
               decMat.real() = c*spasm.mat();
            }
            else
            {
               auto E = nds.find(NonDimensional::Ekman::id())->second->value();

               auto c = -(k1*k1 + k2*k2)/E;
               SparseSM::Chebyshev::LinearMap::I2D1     i2d1(nN, nN, zi, zo);
               decMat.real() = c*i2d1.mat();
            }
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         if(rowId == colId)
         {
            auto laplh = -(k1*k1 + k2*k2);
            if(this->useSplitEquation())
            {
               if(isSplitOperator)
               {
                  SparseSM::Chebyshev::LinearMap::I2Lapl spasm(nN, nN, zi, zo, k1, k2);
                  decMat.real() = laplh*spasm.mat();
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I2Lapl spasm(nN, nN, zi, zo, k1, k2);
                  decMat.real() = laplh*spasm.mat();
               }
            }
            else
            {
               if(k1 == 0 && k2 == 0)
               {
                  SparseSM::Chebyshev::LinearMap::Id spasm(nN, nN, zi, zo, -2, 0);
                  decMat.real() = spasm.mat();
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I4Lapl2 spasm(nN, nN, zi, zo, k1, k2);
                  decMat.real() = laplh*spasm.mat();
               }
            }
         }
         else if(colId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR))
         {
            auto E = nds.find(NonDimensional::Ekman::id())->second->value();

            if(k1 == 0 && k2 == 0)
            {
               auto c = -1./E;
               SparseSM::Chebyshev::LinearMap::I2 i2(nN, nN, zi, zo);
               decMat.real() = c*i2.mat();
            }
            else
            {
               auto c = (k1*k1 + k2*k2)/E;
               SparseSM::Chebyshev::LinearMap::I4D1 i4d1(nN, nN, zi, zo);
               decMat.real() = c*i4d1.mat();
            }
         }
         else if(colId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
         {
            auto Ra = nds.find(NonDimensional::Rayleigh::id())->second->value();
            auto Pr = nds.find(NonDimensional::Prandtl::id())->second->value();

            auto laplh = -(k1*k1 + k2*k2);
            SparseSM::Chebyshev::LinearMap::I4 spasm(nN, nN, zi, zo);
            decMat.real() = -(Ra/Pr)*laplh*spasm.mat();
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && rowId == colId)
      {
         auto Pr = nds.find(NonDimensional::Prandtl::id())->second->value();

         if(k1 == 0 && k2 == 0)
         {
            SparseSM::Chebyshev::LinearMap::Id spasm(nN, nN, zi, zo, -2, 0);
            decMat.real() = (1.0/Pr)*spasm.mat();
         }
         else
         {
            SparseSM::Chebyshev::LinearMap::I2Lapl spasm(nN, nN, zi, zo, k1, k2);
            decMat.real() = (1.0/Pr)*spasm.mat();
         }
      }
      else
      {
         throw std::logic_error("Equations are not setup properly");
      }
   }

   void ModelBackend::timeBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 2);
      int k1 = eigs.at(0);
      int k2 = eigs.at(1);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, k1)(0);

      auto zi = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto zo = nds.find(NonDimensional::Upper1d::id())->second->value();

      if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
      {
         if(k1 == 0 && k2 == 0)
         {
            SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, zi, zo);
            decMat.real() = spasm.mat();
         }
         else
         {
            auto laplh = -(k1*k1 + k2*k2);
            SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, zi, zo);
            decMat.real() = laplh*spasm.mat();
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         if(this->useSplitEquation())
         {
            SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, zi, zo);
            decMat.real() = spasm.mat();
         }
         else
         {
            if(k1 == 0 && k2 == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, zi, zo);
               decMat.real() = spasm.mat();
            }
            else
            {
               auto laplh = -(k1*k1 + k2*k2);
               SparseSM::Chebyshev::LinearMap::I4Lapl spasm(nN, nN, zi, zo, k1, k2);
               decMat.real() = laplh*spasm.mat();
            }
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
      {
         SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, zi, zo);
         decMat.real() = spasm.mat();
      }
   }

   void ModelBackend::splitBoundaryValueBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 2);
      int k1 = eigs.at(0);
      int k2 = eigs.at(1);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, k1)(0);

      if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         decMat.real().resize(nN, 2);
         decMat.imag().resize(nN, 2);

         Eigen::Triplet<MHDFloat> valTop = {0, 0, 1.0};
         Eigen::Triplet<MHDFloat> valBottom = {1, 1, 1.0};
         std::vector<Eigen::Triplet<MHDFloat> > triplets = {valTop, valBottom};
         decMat.real().setFromTriplets(triplets.begin(), triplets.end());
         decMat.imag().setFromTriplets(triplets.begin(), triplets.end());
      }
   }

   void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 2);
      int k1 = eigs.at(0);
      int k2 = eigs.at(1);

      // Time operator
      if(opId == ModelOperator::Time::id())
      {
         bool needStencil = (this->useGalerkin() && bcType == ModelOperatorBoundary::SolverNoTau::id());
         bool needTau = bcType == ModelOperatorBoundary::SolverHasBc::id();

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            this->timeBlock(rModelMatrix, *pRowId, matIdx, res, eigs, nds);

            // Apply boundary condition
            if(needStencil)
            {
               this->applyGalerkinStencil(rModelMatrix.real(), *pRowId, *pRowId, k1, k2, res, bcs, nds);
            }
            else if(needTau)
            {
               this->applyTau(rModelMatrix.real(), *pRowId, *pRowId, k1, k2, res, bcs, nds, false);
            }
         }
      }
      // Linear operator
      else if(opId == ModelOperator::ImplicitLinear::id() || opId == ModelOperator::SplitImplicitLinear::id())
      {
         bool isSplit = (opId == ModelOperator::SplitImplicitLinear::id());
         bool needStencil = (this->useGalerkin() && bcType == ModelOperatorBoundary::SolverNoTau::id());
         bool needTau = bcType == ModelOperatorBoundary::SolverHasBc::id();

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               this->implicitBlock(rModelMatrix, *pRowId, *pColId, matIdx, res, eigs, nds, isSplit);

               // Apply boundary condition
               if(needStencil)
               {
                  this->applyGalerkinStencil(rModelMatrix.real(), *pRowId, *pColId, k1, k2, res, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(rModelMatrix.real(), *pRowId, *pColId, k1, k2, res, bcs, nds, isSplit);
               }
            }
         }
      }
      // Boundary operator
      else if(opId == ModelOperator::Boundary::id() || opId == ModelOperator::SplitBoundary::id())
      {
         bool isSplit = (opId == ModelOperator::SplitBoundary::id());
         bool needStencil = this->useGalerkin();
         bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id());

         auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, k1)(0);

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               rModelMatrix.real().resize(nN, nN);

               // Apply boundary condition
               if(needStencil)
               {
                  this->applyGalerkinStencil(rModelMatrix.real(), *pRowId, *pColId, k1, k2, res, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(rModelMatrix.real(), *pRowId, *pColId, k1, k2, res, bcs, nds, isSplit);
               }
            }
         }
      }
      // Split equation boundary value
      else if(opId == ModelOperator::SplitBoundaryValue::id())
      {
         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            this->splitBoundaryValueBlock(rModelMatrix, *pRowId, matIdx, res, eigs, nds);
         }
      }
      else
      {
         throw std::logic_error("Requested operator type is not implemented");
      }
   }

   void ModelBackend::galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 2);
      int k1 = eigs.at(0);
      int k2 = eigs.at(1);

      this->stencil(mat, fieldId, k1, k2, res, makeSquare, bcs, nds);
   }

   void ModelBackend::explicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const std::size_t opId,  const SpectralFieldId colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 2);
      int k1 = eigs.at(0);
      int k2 = eigs.at(1);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, k1)(0);

      auto zi = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto zo = nds.find(NonDimensional::Upper1d::id())->second->value();

      int s = this->nBc(rowId);
      // Explicit linear operator
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         throw std::logic_error("There are no explicit linear operators");
      }
      // Explicit nonlinear operator
      else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && rowId == colId)
         {
            SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, zi, zo);
            decMat.real() = spasm.mat();
         }
         else
         {
            throw std::logic_error("There are no explicit nonlinear operators");
         }
      }
      // Explicit nextstep operator
      else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         throw std::logic_error("There are no explicit nextstep operators");
      }

      if(this->useGalerkin())
      {
         SparseSM::Chebyshev::LinearMap::Id qId(nN-s, nN, zi, zo, 0, s);
         decMat.real() = qId.mat()*decMat.real();
      }
   }

} // Explicit
} // RRBC
} // Plane
} // Boussinesq
} // Model
} // QuICC
