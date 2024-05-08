/**
 * @file Momentum.hpp
 * @brief Implementation of the vector momentum equation for rotating Rayleigh-Benard convection in a plane layer (toroidal/poloidal formulation)
 */

#ifndef QUICC_EQUATIONS_BOUSSINESQ_PLANE_RRBC_MOMENTUM_HPP
#define QUICC_EQUATIONS_BOUSSINESQ_PLANE_RRBC_MOMENTUM_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

   /**
    * @brief Implementation of the vector momentum equation for rotating Rayleigh-Benard convection in a plane layer (toroidal/poloidal formulation)
    */
   class Momentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams    Shared equation parameters
          */
         Momentum(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Momentum() = default;

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false) override;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements() override;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() override;

         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents() override;

      private:
   };

} // namespace RRBC
} // namespace Plane
} // namespace Boussinesq
} // namespace Equations
} // namepsace QuICC

#endif // QUICC_EQUATIONS_BOUSSINESQ_PLANE_RRBC_MOMENTUM_HPP
