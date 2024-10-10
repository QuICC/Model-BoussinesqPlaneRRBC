/**
 * @file IRRBCBackend.hpp
 * @brief Base model backend for RRBC model
 */

#ifndef QUICC_MODEL_BOUSSINESQ_PLANE_RRBC_IRRBCBACKEND_HPP
#define QUICC_MODEL_BOUSSINESQ_PLANE_RRBC_IRRBCBACKEND_HPP

// System includes
//
#include <map>
#include <memory>
#include <string>
#include <vector>

// Project includes
//
#include "QuICC/Model/IPlaneModelBackend.hpp"
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

/**
 * @brief Base model backed for RRBC model
 */
class IRRBCBackend : public IPlaneModelBackend
{
public:
   /**
    * @brief Constructor
    */
   IRRBCBackend() = default;

   /**
    * @brief Destructor
    */
   virtual ~IRRBCBackend() = default;

   /**
    * @brief Get vector of names for the physical fields
    */
   virtual std::vector<std::string> fieldNames() const override;

   /**
    * @brief Get vector of names for the nondimensional parameters
    */
   virtual std::vector<std::string> paramNames() const override;

   /**
    * @brief Get vector of bools about periodic box
    */
   virtual std::vector<bool> isPeriodicBox() const override;

   /**
    * @brief Get automatically computed parameters based on input parameters
    *
    * @param cfg  Input parameters
    */
   virtual std::map<std::string, MHDFloat> automaticParameters(
      const std::map<std::string, MHDFloat>& cfg) const override;

protected:
   /**
    * @brief Number of boundary conditions
    *
    * @fId  Field ID
    */
   int nBc(const SpectralFieldId& fId) const override;

   /**
    * @brief Apply tau line for boundary condition
    *
    * @param mat     Input/Output matrix to apply tau line to
    * @param rowId   ID of field of equation
    * @param colId   ID of field
    * @param k1      First wave number
    * @param opts    Additional options
    * @param res     Resolution object
    * @param bcs     Boundary conditions
    * @param nds     Nondimensional parameters
    * @param isSplitOperator  Is second operator of split 4th order system?
    */
   void applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const int k1,
      std::shared_ptr<details::BlockOptions> opts, const Resolution& res,
      const BcMap& bcs, const NonDimensional::NdMap& nds,
      const bool isSplitOperator) const override;

   /**
    * @brief Boundary condition stencil
    *
    * @param mat        Input/Output matrix to store galerkin stencil
    * @param fID        Field ID
    * @param k1         First wave number
    * @param res        Resolution object
    * @param makeSquare Truncate operator to make square
    * @param bcs        Boundary conditions
    * @param nds        Nondimensional parameters
    */
   virtual void stencil(SparseMatrix& mat, const SpectralFieldId& fId,
      const int k1, const Resolution& res, const bool makeSquare,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const;

   /**
    * @brief Apply galerkin stencil for boundary condition
    *
    * @param mat     Input/Output matrix to apply stencil to
    * @param rowId   ID of field of equation
    * @param colId   ID of field
    * @param k1r     Row space first wave number
    * @param k1c     Column space first wave number
    * @param opts    Additional options
    * @param res     Resolution object
    * @param bcs     Boundary conditions
    * @param nds     Nondimensional parameters
    */
   void applyGalerkinStencil(SparseMatrix& decMat, const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const int k1r, const int k1c,
      std::shared_ptr<details::BlockOptions> opts, const Resolution& res,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

private:
};

namespace implDetails {

/**
 * @brief Specific options for current model
 */
struct BlockOptionsImpl : public details::BlockOptions
{
   /**
    * @brief default ctor
    */
   BlockOptionsImpl() = default;

   /**
    * @brief default dtor
    */
   virtual ~BlockOptionsImpl() = default;

   /// Lower boundary
   Internal::MHDFloat zi;
   /// Upper boundary
   Internal::MHDFloat zo;
   /// First wave number
   MHDFloat k1;
   /// Second wave number
   MHDFloat k2;
   /// Use truncated quasi-inverse?
   bool truncateQI;
   /// Boundary condition
   std::size_t bcId;
   /// Split operator for influence matrix?
   bool isSplitOperator;
   /// Use split equation for influence matrix?
   bool useSplitEquation;
};
} // namespace implDetails

} // namespace RRBC
} // namespace Plane
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_PLANE_RRBC_IRRBCBACKEND_HPP
