/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating Rayleigh-Benard in a plane
 * layer (toroidal/poloidal formulation) with anisotropic rescaling model
 */

#ifndef QUICC_MODEL_BOUSSINESQ_PLANE_RRBC_EXPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_PLANE_RRBC_EXPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/SpatialScheme/3D/TFF.hpp"
#include "QuICC/Model/PyModelBackend.hpp"
#include "Model/Boussinesq/Plane/RRBC/Explicit/ModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

namespace Explicit {

/**
 * @brief Implementation of the Boussinesq rotating Rayleigh-Benard in a plane
 * layer (toroidal/poloidal formulation) with anisotropic rescaling model
 */
template <typename TBuilder> class PhysicalModel : public TBuilder
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

   /**
    * @brief Initialize specialized backend
    */
   void init() final;

protected:
private:
};

template <typename TBuilder> void PhysicalModel<TBuilder>::init()
{
   TBuilder::init();
#ifdef QUICC_MODEL_BOUSSINESQPLANERRBC_EXPLICIT_BACKEND_CPP

   this->mpBackend = std::make_shared<ModelBackend>();
#else
   std::string pyModule = "boussinesq.plane.rrbc.explicit.physical_model";
   std::string pyClass = "PhysicalModel";

   this->mpBackend =
      std::make_shared<PyModelBackend>(pyModule, pyClass);
#endif
}

} // namespace Explicit
} // namespace RRBC
} // namespace Plane
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC


#endif // QUICC_MODEL_BOUSSINESQ_PLANE_RRBC_EXPLICIT_PHYSICALMODEL_HPP
