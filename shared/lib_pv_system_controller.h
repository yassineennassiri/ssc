#ifndef _LIB_PV_SYSTEM_CONTROLLER_H_
#define _LIB_PV_SYSTEM_CONTROLLER_H_

#include <memory>

#include "lib_pv_io.h"

class Subarray;
class MPPTController;
class PVDCController;
class PVACController;
class PVSystem;

/**
* \class PVSystemController
*
* This class oversees the running of the PVSystem and other potential technologies running independently of the PVSystem (hybrid system)
* As an example, the PVSystemController would run the PVSystem for one time step, run a WindSystem for one timestep, and then combine the output
* It is broadly intended as a supervisory controller for complex system configurations with multiple technologies
*
*/
class PVSystemController
{
public:
	/// Construct a PVSystemController from a PVIOManager
	PVSystemController(std::shared_ptr<PVIOManager> pvIOManager);

	/// Run the System for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a shared pointer
	std::shared_ptr<PVIOManager> m_pvIOManager;

	/// The PVSystemController exclusively manages the PVSystem
	std::unique_ptr<PVSystem> m_pvSystem;
};

/**
* \class PVSystem
*
* A PVSystem contains all of the components, models, and controllers needed to operate a PV plus storage system configured 
* on the AC or DC side of the array.
*
*/
class PVSystem
{
public:

	/// Construct a PVSystem from a PVIOManager
	PVSystem(std::shared_ptr<PVIOManager> pvIOManager);

	/// Simulation the PVSystem for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a shared pointer
	std::shared_ptr<PVIOManager> m_pvIOManager;

	/// The PVSystem uniquely manages the DC and AC controllers
	std::unique_ptr<PVDCController> m_pvDCController;
	std::unique_ptr<PVACController> m_pvACController;
};

/**
* \class PVDCController
*
* A PVDCController contains the components and models required to operate all of the DC components of a PV array
* including all subarray's, and DC-connected batteries
*
*/
class PVDCController
{
public: 

	/// Construct a PVDCController from a PVIOManager
	PVDCController(std::shared_ptr<PVIOManager> pvIOManager);

	//const bool RunSingleStep();
private:

	/// The PVIOManager is a shared pointer
	std::shared_ptr<PVIOManager> m_pvIOManager;

	/// The PVDCController uniquely manages MPPTControllers in the system
	std::vector<std::unique_ptr<MPPTController>> m_MPPTControllers;

};

/**
* \class PVACController
*
* A PVACController contains the components and models required to operate all of the AC components of a PV array
* including all subarray's, and an AC-connected battery
*
*/
class PVACController
{
public:
	/// Construct a PVDAController from a PVIOManager
	PVACController(std::shared_ptr<PVIOManager> pvIOManager);

	//const bool RunSingleStep();

private:
	/// The PVIOManager is a shared pointer
	std::shared_ptr<PVIOManager> m_pvIOManager;

};

/**
* \class MPPTController
*
* A MPPTController contains information around MPPT ports and their operation limits, and information about
* which Subbarray's are controlled
*
*/
class MPPTController
{
public:
	/// Construct a MPPTController from a PVIOManager
	MPPTController(std::shared_ptr<PVIOManager> pvIOManager);

	/// Simulation the MPPTController for one time step
	const bool RunSingleStep();

private:
	/// The PVIOManager is a shared pointer
	std::shared_ptr<PVIOManager> m_pvIOManager;

	/// The MPPTController uniquely manages the Subarray's under its control
	std::vector<std::unique_ptr<Subarray>> m_Subarrays;

	/// The MPPTController_IO is part of the PVIOManager, but broken out for convenience
	MPPTController_IO * m_MPPTIO;
};

/**
* \class Subarray
*
* A Subarray contains information relevent to one subarray, including which irradiance model should be used,
* which module model should be used, what the module level properties are, and whether there is a battery 
* connected the subarray on the DC side.
*
*/
class Subarray
{
public:
	/// Construct a Subarray from a PVIOManager
	Subarray(std::shared_ptr<PVIOManager> pvIOManager);

	/// Simulate a Subarray for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a shared pointer
	std::shared_ptr<PVIOManager> m_pvIOManager;
	
};

#endif
