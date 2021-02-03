/** 
 * @file Model.cpp
 * @brief Source of the rotating Boussinesq Rayleigh-Benard convection in a plane layer (toroidal/poloidal formulation) model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/RRBC/PhysicalModel.hpp )

// Project includes
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/RRBC/Transport.hpp )
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/RRBC/Momentum.hpp )
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Io/Variable/Cartesian1DScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/Cartesian1DTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/Cartesian1DNusseltDZWriter.hpp"
#include "QuICC/Io/Variable/VisualizationFileWriter.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Generator/States/RandomScalarState.hpp"
#include "QuICC/Generator/States/RandomVectorState.hpp"
#include "QuICC/Generator/States/CartesianExactScalarState.hpp"
#include "QuICC/Generator/States/CartesianExactVectorState.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldTrivialVisualizer.hpp"
#include "QuICC/Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "QuICC/Model/PhysicalModelBase.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

   const std::string PhysicalModel::PYMODULE = "boussinesq_rrbcplane";

   const std::string PhysicalModel::PYCLASS = "BoussinesqRRBCPlane";

   void PhysicalModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addScalarEquation<Equations::Boussinesq::Plane::RRBC::Transport>();
      
      // Add Navier-Stokes equation
      spSim->addVectorEquation<Equations::Boussinesq::Plane::RRBC::Momentum>();
   }

   void PhysicalModel::addStates(SharedStateGenerator spGen)
   {
      // Generate "exact" solutions (trigonometric or monomial)
      if(false)
      {
         // Shared pointer to equation
         Equations::SharedCartesianExactScalarState spScalar;
         Equations::SharedCartesianExactVectorState spVector;

         // Add vector exact initial state generator
         spVector = spGen->addVectorEquation<Equations::CartesianExactVectorState>();
         spVector->setIdentity(PhysicalNames::Velocity::id());
         spVector->setStateType(Equations::CartesianExactStateIds::TORPOLTFF);

         // Add scalar exact initial state generator
         spScalar = spGen->addScalarEquation<Equations::CartesianExactScalarState>();
         spScalar->setIdentity(PhysicalNames::Temperature::id());
         spScalar->setStateType(Equations::CartesianExactStateIds::PLANFORMSQUARES);
         spScalar->setModeOptions(1e0, 10.0, 1e0, 10.0);

      // Generate random spectrum
      } else
      {
         // Shared pointer to random initial state equation
         Equations::SharedRandomScalarState spScalar;
         Equations::SharedRandomVectorState spVector;

         // Add scalar random initial state generator 
         spVector = spGen->addVectorEquation<Equations::RandomVectorState>();
         spVector->setIdentity(PhysicalNames::Velocity::id());
         spVector->setSpectrum(FieldComponents::Spectral::TOR, -1e-6, 1e-6, 1e4, 1e4, 1e4);
         spVector->setSpectrum(FieldComponents::Spectral::POL, -1e-6, 1e-6, 1e4, 1e4, 1e4);

         // Add scalar random initial state generator
         spScalar = spGen->addScalarEquation<Equations::RandomScalarState>();
         spScalar->setIdentity(PhysicalNames::Temperature::id());
         spScalar->setSpectrum(-1e-6, 1e-6, 1e4, 1e4, 1e4);
      }

      // Add output file
      Io::Variable::SharedStateFileWriter spOut(new Io::Variable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      spOut->expect(PhysicalNames::Velocity::id());
      spOut->expect(PhysicalNames::Temperature::id());
      spGen->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedScalarFieldTrivialVisualizer spSTrivial;
      Equations::SharedVectorFieldVisualizer spVector;

      // Add temperature field visualization
      spScalar = spVis->addScalarEquation<Equations::ScalarFieldVisualizer>();
      spScalar->setFields(true, false);
      spScalar->setIdentity(PhysicalNames::Temperature::id());

      // Add mean temperature field visualization
      spSTrivial = spVis->addScalarEquation<Equations::ScalarFieldTrivialVisualizer>();
      spSTrivial->setFields(true, false);
      spSTrivial->setIdentity(PhysicalNames::MeanTemperature::id());

      // Add fluctuating temperature field visualization
      spSTrivial = spVis->addScalarEquation<Equations::ScalarFieldTrivialVisualizer>();
      spSTrivial->setFields(true, false);
      spSTrivial->setIdentity(PhysicalNames::FluctTemperature::id());

      // Add velocity fields visualization
      spVector = spVis->addVectorEquation<Equations::VectorFieldVisualizer>();
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::Velocity::id());

      // Add output file
      Io::Variable::SharedVisualizationFileWriter spOut(new Io::Variable::VisualizationFileWriter(SchemeType::type()));
      spOut->expect(PhysicalNames::Temperature::id());
      spOut->expect(PhysicalNames::MeanTemperature::id());
      spOut->expect(PhysicalNames::FluctTemperature::id());
      spOut->expect(PhysicalNames::Velocity::id());
      spVis->addHdf5OutputFile(spOut);
   }

   void PhysicalModel::setVisualizationState(SharedVisualizationGenerator spVis)
   {
      // Create and add initial state file to IO
      Io::Variable::SharedStateFileReader spIn(new Io::Variable::StateFileReader("4Visu", SchemeType::type(), SchemeType::isRegular()));

      // Set expected fields
      spIn->expect(PhysicalNames::Temperature::id());
      spIn->expect(PhysicalNames::Velocity::id());

      // Set simulation state
      spVis->setInitialState(spIn);
   }

   void PhysicalModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create temperature energy writer
      Io::Variable::SharedCartesian1DScalarEnergyWriter spTemp(new Io::Variable::Cartesian1DScalarEnergyWriter("temperature", SchemeType::type()));
      spTemp->expect(PhysicalNames::Temperature::id());
      spSim->addAsciiOutputFile(spTemp);

      // Create kinetic energy writer
      Io::Variable::SharedCartesian1DTorPolEnergyWriter spKinetic(new Io::Variable::Cartesian1DTorPolEnergyWriter("kinetic", SchemeType::type()));
      spKinetic->expect(PhysicalNames::Velocity::id());
      spSim->addAsciiOutputFile(spKinetic);

      // Create nusselt number writer
      Io::Variable::SharedCartesian1DNusseltDZWriter spNusselt(new Io::Variable::Cartesian1DNusseltDZWriter(SchemeType::type()));
      spNusselt->expect(PhysicalNames::Temperature::id());
      spSim->addAsciiOutputFile(spNusselt);
   }

   void PhysicalModel::addHdf5OutputFiles(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<std::size_t>::const_iterator  it;
      std::vector<std::size_t> ids = PhysicalModelBase::fieldIds();

      // Create and add state file to IO
      Io::Variable::SharedStateFileWriter spState(new Io::Variable::StateFileWriter(SchemeType::type(), SchemeType::isRegular()));
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spState->expect(*it);
      }

      spSim->addHdf5OutputFile(spState);
   }

   void PhysicalModel::addStatsOutputFiles(SharedSimulation spSim)
   {
   }

   void PhysicalModel::setInitialState(SharedSimulation spSim)
   {
      // Field IDs iterator
      std::vector<std::size_t>::const_iterator  it;
      std::vector<std::size_t> ids = PhysicalModelBase::fieldIds();

      // Create and add initial state file to IO
      Io::Variable::SharedStateFileReader spInit(new Io::Variable::StateFileReader("_initial", SchemeType::type(), SchemeType::isRegular()));

      // Set expected field names
      for(it = ids.begin(); it != ids.end(); ++it)
      {
         spInit->expect(*it);
      }

      // Set simulation state
      spSim->setInitialState(spInit);
   }

}
}
}
}
}
