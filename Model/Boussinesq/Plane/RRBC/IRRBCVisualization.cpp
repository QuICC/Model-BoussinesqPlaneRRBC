/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq Rotating Rayleigh-Benard convection in a
 * plane layer (toroidal/poloidal formulation) model
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Plane/RRBC/IRRBCVisualization.hpp"
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

VectorFormulation::Id IRRBCVisualization::SchemeFormulation()
{
   return VectorFormulation::TORPOL;
}

std::string IRRBCVisualization::version() const
{
   return std::string(gitHash);
}

void IRRBCVisualization::addVisualizers(SharedVisualizationGenerator spVis)
{
   // Shared pointer to basic field visualizer
   Equations::SharedScalarFieldVisualizer spScalar;
   Equations::SharedVectorFieldVisualizer spVector;

   // Add temperature field visualization
   spScalar =
      spVis->addEquation<Equations::ScalarFieldVisualizer>(this->spBackend());
   spScalar->setFields(true, false);
   spScalar->setIdentity(PhysicalNames::Temperature::id());

   // Add velocity fields visualization
   spVector =
      spVis->addEquation<Equations::VectorFieldVisualizer>(this->spBackend());
   spVector->setFields(true, false, true);
   spVector->setIdentity(PhysicalNames::Velocity::id());

   // Add output file
   auto spOut = std::make_shared<Io::Variable::VisualizationFileWriter>(
      spVis->ss().tag());
   spOut->expect(PhysicalNames::Temperature::id());
   spOut->expect(PhysicalNames::Velocity::id());
   spVis->addHdf5OutputFile(spOut);
}

} // namespace RRBC
} // namespace Plane
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
