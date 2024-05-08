/** 
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating Rayleigh-Benard in a plane layer (toroidal/poloidal formulation) with anisotropic rescaling model
 */

#ifndef QUICC_MODEL_BOUSSINESQ_PLANE_RRBC_EXPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_PLANE_RRBC_EXPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "Model/Boussinesq/Plane/RRBC/IRRBCModel.hpp"
#include "QuICC/SpatialScheme/3D/TFF.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

namespace Explicit {

   /**
    * @brief Implementation of the Boussinesq rotating Rayleigh-Benard in a plane layer (toroidal/poloidal formulation) with anisotropic rescaling model
    */
   class PhysicalModel: public IRRBCModel
   {
      public:
         /// Typedef for the spatial scheme used
         typedef SpatialScheme::TFF SchemeType;

         /**
          * @brief Constructor
          */
         PhysicalModel() = default;

         /**
          * @brief Destructor
          */
         virtual ~PhysicalModel() = default;

         /// Python model script module name
         virtual std::string PYMODULE() override;

         /**
          * @brief Initialize specialized backend
          */
         void init() final;

      protected:

      private:
   };

} // namespace Explicit
} // namespace RRBC
} // namespace Plane
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC


#endif // QUICC_MODEL_BOUSSINESQ_PLANE_RRBC_EXPLICIT_PHYSICALMODEL_HPP
