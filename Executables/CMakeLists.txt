include(ModelFunctions)

######################
# Explicit executables
set(Type "Explicit")
set(ModelId "${QUICC_CURRENT_MODEL_DIRNAME}/${Type}")
string(TOLOWER "quicc_${QUICC_CURRENT_MODEL_DIRNAME}" _model_lib)
string(TOLOWER "${_model_lib}_${Type}" modLib)
set(PhysicalModel "Model/Boussinesq/Plane/RRBC/Explicit/PhysicalModel.hpp")
set(ClassPhysicalModel "QuICC::Model::Boussinesq::Plane::RRBC::Explicit::PhysicalModel")

# Config
set(Exe "Config")
set(Factory "QuICC/Model/ConfigFactory.hpp")
set(ClassFactory "QuICC::ConfigFactory")
set(Builder "Model/Boussinesq/Plane/RRBC/IRRBCModel.hpp")
set(ClassBuilder "QuICC::Model::Boussinesq::Plane::RRBC::IRRBCModel")
set(source_ "${CMAKE_CURRENT_BINARY_DIR}/${Type}/Run${Exe}.cpp")
configure_file(RunExe.cpp.in ${source_} @ONLY)

quicc_setup_exe(${ModelId}
  POSTFIX "${Exe}"
  EXESOURCE "${source_}"
  MODELLIB "${modLib}")

# State
set(Exe "State")
set(Factory "QuICC/Model/StateGeneratorFactory.hpp")
set(ClassFactory "QuICC::StateGeneratorFactory")
set(Builder "Model/Boussinesq/Plane/RRBC/IRRBCState.hpp")
set(ClassBuilder "QuICC::Model::Boussinesq::Plane::RRBC::IRRBCState")
set(source_ "${CMAKE_CURRENT_BINARY_DIR}/${Type}/Run${Exe}.cpp")
configure_file(RunExe.cpp.in ${source_} @ONLY)

quicc_setup_exe(${ModelId}
  POSTFIX "${Exe}"
  EXESOURCE "${source_}"
  MODELLIB "${modLib}")

# Model
set(Exe "Model")
set(Factory "QuICC/Model/ModelFactory.hpp")
set(ClassFactory "QuICC::ModelFactory")
set(Builder "Model/Boussinesq/Plane/RRBC/IRRBCModel.hpp")
set(ClassBuilder "QuICC::Model::Boussinesq::Plane::RRBC::IRRBCModel")
set(source_ "${CMAKE_CURRENT_BINARY_DIR}/${Type}/Run${Exe}.cpp")
configure_file(RunExe.cpp.in ${source_} @ONLY)

quicc_setup_exe(${ModelId}
  POSTFIX "${Exe}"
  EXESOURCE "${source_}"
  MODELLIB "${modLib}")

# Visu
set(Exe "Visu")
set(Factory "QuICC/Model/VisualizationGeneratorFactory.hpp")
set(ClassFactory "QuICC::VisualizationGeneratorFactory")
set(Builder "Model/Boussinesq/Plane/RRBC/IRRBCVisualization.hpp")
set(ClassBuilder "QuICC::Model::Boussinesq::Plane::RRBC::IRRBCVisualization")
set(source_ "${CMAKE_CURRENT_BINARY_DIR}/${Type}/Run${Exe}.cpp")
configure_file(RunExe.cpp.in ${source_} @ONLY)

quicc_setup_exe(${ModelId}
  POSTFIX "${Exe}"
  EXESOURCE "${source_}"
  MODELLIB "${modLib}")


######################
# Stability executables
if(TARGET QuICC::Stability)
  set(Type "Explicit")
  set(ModelId "${QUICC_CURRENT_MODEL_DIRNAME}/${Type}")
  string(TOLOWER "quicc_${QUICC_CURRENT_MODEL_DIRNAME}" _model_lib)
  string(TOLOWER "${_model_lib}_${Type}" modLib)
  set(PhysicalModel "Model/Boussinesq/Plane/RRBC/Explicit/PhysicalModel.hpp")
  set(ClassPhysicalModel "QuICC::Model::Boussinesq::Plane::RRBC::Explicit::PhysicalModel")

  # Config
  set(Exe "StabilityConfig")
  set(Factory "Stability/ConfigFactory.hpp")
  set(ClassFactory "QuICC::ConfigFactory")
  set(Builder "Model/Boussinesq/Plane/RRBC/IRRBCModel.hpp")
  set(ClassBuilder "QuICC::Model::Boussinesq::Plane::RRBC::IRRBCModel")
  set(source_ "${CMAKE_CURRENT_BINARY_DIR}/${Type}/Run${Exe}.cpp")
  configure_file(RunExe.cpp.in ${source_} @ONLY)

  quicc_setup_exe(${ModelId}
    POSTFIX "${Exe}"
    EXESOURCE "${source_}"
    MODELLIB "${modLib}"
    EXTRALIBS QuICC::Stability)

  # Stability
  set(Exe "Stability")
  set(Factory "Stability/LinearStabilityFactory.hpp")
  set(ClassFactory "QuICC::LinearStabilityFactory")
  set(Builder "Model/Boussinesq/Plane/RRBC/IRRBCModel.hpp")
  set(ClassBuilder "QuICC::Model::Boussinesq::Plane::RRBC::IRRBCModel")
  set(source_ "${CMAKE_CURRENT_BINARY_DIR}/${Type}/Run${Exe}.cpp")
  configure_file(RunExe.cpp.in ${source_} @ONLY)

  quicc_setup_exe(${ModelId}
    POSTFIX "${Exe}"
    EXESOURCE "${source_}"
    MODELLIB "${modLib}"
    EXTRALIBS QuICC::Stability)
endif()
