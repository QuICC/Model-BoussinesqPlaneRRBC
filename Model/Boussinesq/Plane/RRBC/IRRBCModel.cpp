/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq Rotating Rayleigh-Benard convection in a
 * plane layer (toroidal/poloidal formulation) model
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Plane/RRBC/IRRBCModel.hpp"
#include "Model/Boussinesq/Plane/RRBC/Momentum.hpp"
#include "Model/Boussinesq/Plane/RRBC/Transport.hpp"
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

VectorFormulation::Id IRRBCModel::SchemeFormulation()
{
   return VectorFormulation::TORPOL;
}

std::string IRRBCModel::version() const
{
   return "BoussinesqPlaneRRBC:" + std::string(gitHash);
}

void IRRBCModel::addEquations(SharedSimulation spSim)
{
   // Add transport equation
   spSim->addEquation<Equations::Boussinesq::Plane::RRBC::Transport>(
      this->spBackend());

   // Add Navier-Stokes equation
   spSim->addEquation<Equations::Boussinesq::Plane::RRBC::Momentum>(
      this->spBackend());
}

std::map<std::string, std::map<std::string, int>> IRRBCModel::configTags() const
{
   std::map<std::string, int> onOff;
   onOff.emplace("enable", 1);

   std::map<std::string, int> offOn;
   onOff.emplace("enable", 0);

   std::map<std::string, std::map<std::string, int>> tags;
   // kinetic
   tags.emplace("kinetic_energy", onOff);
   // temperature
   tags.emplace("temperature_energy", onOff);
   tags.emplace("temperature_nusselt", offOn);

   return tags;
}

void IRRBCModel::addAsciiOutputFiles(SharedSimulation spSim)
{
   // Create temperature energy writer
   this->enableAsciiFile<Io::Variable::Cartesian1DScalarEnergyWriter>(
      "temperature_energy", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create kinetic energy writer
   this->enableAsciiFile<Io::Variable::Cartesian1DTorPolEnergyWriter>(
      "kinetic_energy", "kinetic", PhysicalNames::Velocity::id(), spSim);
}

} // namespace RRBC
} // namespace Plane
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
