/**
 * @file Transport.cpp
 * @brief Source of the implementation of the transport equation for rotating
 * Rayleigh-Benard convection in a plane layer (toroidal/poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Plane/RRBC/Transport.hpp"
#include "Model/Boussinesq/Plane/RRBC/TransportKernel.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Transform/Path/ScalarNl.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

Transport::Transport(SharedEquationParameters spEqParams,
   SpatialScheme::SharedCISpatialScheme spScheme,
   std::shared_ptr<Model::IModelBackend> spBackend) :
    IScalarEquation(spEqParams, spScheme, spBackend)
{
   // Set the variable requirements
   this->setRequirements();
}

void Transport::setCoupling()
{
   auto features = defaultCouplingFeature();
   features.at(CouplingFeature::Nonlinear) = true;

   this->defineCoupling(FieldComponents::Spectral::SCALAR,
      CouplingInformation::PROGNOSTIC, 0, features);
}

void Transport::setNLComponents()
{
   this->addNLComponent(FieldComponents::Spectral::SCALAR,
      Transform::Path::ScalarNl::id());
}

void Transport::initNLKernel(const bool force)
{
   // Initialize if empty or forced
   if (force || !this->mspNLKernel)
   {
      // Initialize the physical kernel
      auto spNLKernel = std::make_shared<Physical::Kernel::TransportKernel>();
      spNLKernel->setScalar(this->name(), this->spUnknown());
      spNLKernel->setVector(PhysicalNames::Velocity::id(),
         this->spVector(PhysicalNames::Velocity::id()));
      spNLKernel->init(1.0);
      this->mspNLKernel = spNLKernel;
   }
}

void Transport::setRequirements()
{
   // Set temperatur as equation unknown
   this->setName(PhysicalNames::Temperature::id());

   // Set solver timing
   this->setSolveTiming(SolveTiming::Prognostic::id());

   // Forward transform generates nonlinear RHS
   this->setForwardPathsType(FWD_IS_NONLINEAR);

   // Get reference to spatial scheme
   const auto& ss = this->ss();

   // Add temperature to requirements: is scalar?, need spectral?, need
   // physical?, need diff?
   auto& tempReq =
      this->mRequirements.addField(PhysicalNames::Temperature::id(),
         FieldRequirement(true, ss.spectral(), ss.physical()));
   tempReq.enableSpectral();
   tempReq.enableGradient();

   // Add X velocity to requirements: is scalar?, need spectral?, need
   // physical?, need diff?
   auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(),
      FieldRequirement(false, ss.spectral(), ss.physical()));
   velReq.enableSpectral();
   velReq.enablePhysical();
}

} // namespace RRBC
} // namespace Plane
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC
