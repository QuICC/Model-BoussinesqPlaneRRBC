/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector momentum equation for rotating Rayleigh-Benard convection in a plane layer (toroidal/poloidal formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Model/Boussinesq/Plane/RRBC/Momentum.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

   Momentum::Momentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Momentum::~Momentum()
   {
   }

   void Momentum::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, 0, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, 0, true, false);
   }

   void Momentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, 0);

      this->addNLComponent(FieldComponents::Spectral::POL, 0);
   }

   void Momentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      ///
      /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
      ///
      switch(compId)
      {
         case(FieldComponents::Physical::X):
            Physical::Cross<FieldComponents::Physical::Y,FieldComponents::Physical::Z>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::Y):
            Physical::Cross<FieldComponents::Physical::Z,FieldComponents::Physical::X>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::Z):
            Physical::Cross<FieldComponents::Physical::X,FieldComponents::Physical::Y>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            break;
         default:
            assert(false);
            break;
      }
   }

   void Momentum::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need grad?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, true));
   }

}
}
}
}
}
