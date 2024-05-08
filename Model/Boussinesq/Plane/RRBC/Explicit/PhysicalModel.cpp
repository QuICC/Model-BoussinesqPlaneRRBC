/** 
 * @file Model.cpp
 * @brief Source of the rotating Boussinesq Rayleigh-Benard convection in a plane layer (toroidal/poloidal formulation) model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Plane/RRBC/Explicit/PhysicalModel.hpp"
#include "Model/Boussinesq/Plane/RRBC/Explicit/ModelBackend.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Plane {

namespace RRBC {

namespace Explicit {

   std::string PhysicalModel::PYMODULE()
   {
      return "boussinesq.plane.rrbc.explicit.physical_model";
   }

   void PhysicalModel::init()
   {
#ifdef QUICC_MODEL_BOUSSINESQPLANERRBC_EXPLICIT_BACKEND_CPP
      IPhysicalModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<ModelBackend>();
#else
      IPhysicalPyModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
#endif
   }

}
}
}
}
}
}
