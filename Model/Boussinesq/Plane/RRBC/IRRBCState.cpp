/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq Rotating Rayleigh-Benard convection in a
 * plane layer (toroidal/poloidal formulation) model
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Plane/RRBC/IRRBCState.hpp"
#include "Model/Boussinesq/Plane/RRBC/gitHash.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/Cartesian1DScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/Cartesian1DTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Io/Variable/VisualizationFileWriter.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Generator/States/CartesianExactScalarState.hpp"
#include "QuICC/Generator/States/CartesianExactVectorState.hpp"
#include "QuICC/Generator/States/RandomScalarState.hpp"
#include "QuICC/Generator/States/RandomVectorState.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldTrivialVisualizer.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "QuICC/SpectralKernels/MakeRandom.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

VectorFormulation::Id IRRBCState::SchemeFormulation()
{
   return VectorFormulation::TORPOL;
}

std::string IRRBCState::version() const
{
   return std::string(gitHash);
}

void IRRBCState::addStates(SharedStateGenerator spGen)
{
   // Shared pointer to equation
   Equations::SharedCartesianExactScalarState spScalar;
   Equations::SharedCartesianExactVectorState spVector;

   // Add temperature initial state generator
   spScalar = spGen->addEquation<Equations::CartesianExactScalarState>(
      this->spBackend());
   spScalar->setIdentity(PhysicalNames::Temperature::id());
   switch (0)
   {
   case 0: {
      spScalar->setPhysicalNoise(1e-15);
   }
   break;

   case 1: {
      spScalar->setPhysicalConstant(1.0);
   }
   break;

   case 3: {
      auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(
         spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
      std::vector<MHDFloat> ratios = {1e2, 1e2, 1e2};
      spKernel->setRatio(ratios);
      spKernel->init(-1e-15, 1e-15);
      spVector->setSrcKernel(FieldComponents::Spectral::SCALAR, spKernel);
   }
   break;
   }

   // Add velocity initial state generator
   spVector = spGen->addEquation<Equations::CartesianExactVectorState>(
      this->spBackend());
   spVector->setIdentity(PhysicalNames::Velocity::id());
   switch (2)
   {
   case 0: {
      spScalar->setPhysicalNoise(1e-15);
   }
   break;

   case 1: {
      spScalar->setPhysicalConstant(1.0);
   }
   break;

   case 2: {
      auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(
         spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
      std::vector<MHDFloat> ratios = {1e2, 1e2, 1e2};
      spKernel->setRatio(ratios);
      spKernel->init(-1e-15, 1e-15);
      spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
      spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
   }
   break;

   case 3: {
      auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(
         spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
      std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
      spKernel->setRatio(ratios);
      spKernel->init(-1e-7, 1e-7);
      spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
      spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
   }
   break;

   case 4: {
      // spVector->setStateType(Equations::CartesianExactStateIds::TORPOLTFF);
   }
   break;
   }

   // Add output file
   auto spOut =
      std::make_shared<Io::Variable::StateFileWriter>(spGen->ss().tag(),
         spGen->ss().has(SpatialScheme::Feature::RegularSpectrum));
   spOut->expect(PhysicalNames::Velocity::id());
   spOut->expect(PhysicalNames::Temperature::id());
   spGen->addHdf5OutputFile(spOut);
}

} // namespace RRBC
} // namespace Plane
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
