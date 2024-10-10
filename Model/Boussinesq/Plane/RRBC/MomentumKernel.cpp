/**
 * @file MomentumKernel.cpp
 * @brief Source of physical space kernel for the Momentum equation
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Plane/RRBC/MomentumKernel.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

MomentumKernel::MomentumKernel() : IPhysicalKernel() {}

std::size_t MomentumKernel::name() const
{
   return this->mName;
}

void MomentumKernel::setVelocity(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mName = name;

   this->setField(name, spField);
}

void MomentumKernel::init(const MHDFloat inertia)
{
   // Set scaling constants
   this->mInertia = inertia;
}

void MomentumKernel::compute(Framework::Selector::PhysicalScalarField& rNLComp,
   FieldComponents::Physical::Id id) const
{
   ///
   /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
   ///
   switch (id)
   {
   case (FieldComponents::Physical::X):
      std::visit(
         [&](auto&& v)
         {
            Physical::Cross<FieldComponents::Physical::Y,
               FieldComponents::Physical::Z>::set(rNLComp, v->dom(0).curl(),
               v->dom(0).phys(), this->mInertia);
         },
         this->vector(this->name()));
      break;
   case (FieldComponents::Physical::Y):
      std::visit(
         [&](auto&& v)
         {
            Physical::Cross<FieldComponents::Physical::Z,
               FieldComponents::Physical::X>::set(rNLComp, v->dom(0).curl(),
               v->dom(0).phys(), this->mInertia);
         },
         this->vector(this->name()));
      break;
   case (FieldComponents::Physical::Z):
      std::visit(
         [&](auto&& v)
         {
            Physical::Cross<FieldComponents::Physical::X,
               FieldComponents::Physical::Y>::set(rNLComp, v->dom(0).curl(),
               v->dom(0).phys(), this->mInertia);
         },
         this->vector(this->name()));
      break;
   default:
      assert(false);
      break;
   }
}

} // namespace Kernel
} // namespace Physical
} // namespace QuICC
