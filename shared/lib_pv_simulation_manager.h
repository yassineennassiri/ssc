#ifndef _LIB_PV_SIMULATION_MANAGER_H_
#define _LIB_PV_SIMULATION_MANAGER_H_

#include <memory>

#include "lib_pv_io.h"

class PVLossManager;
class PVSystemController;

/** 
* \class PVSimulationManager
*
*  PVSimulationManager manages PV simulation details, including managing timestep information,
*  and invoking the SystemController.  The manager contains a PVIOManager and PVLossManager
*/
class PVSimulationManager
{
public: 

	/// Construct a PVSimulationManager with a PVIOManager and PVLossManager
	PVSimulationManager(compute_module &cm);

	/// PVSimulationManager is intended to be non-copyable
	PVSimulationManager(const PVSimulationManager&);
	PVSimulationManager& operator=(const PVSimulationManager&);

	/// Run a PV Simulation
	const bool Simulate();

	/// Return a weak pointer to the PVIOManager
	PVIOManager * getPVIOManager() { return m_PVIOManager.get(); }

	/// Return a weak pointer to the PVLossManager
	PVLossManager * getPVLossManager() { return m_PVLossManager.get(); }


private:

	/// Run a single step of the simulation
	const bool RunSingleStep();

	// These objects are managed exclusively by PVSimulationManager
	std::unique_ptr<PVIOManager> m_PVIOManager;			/// An object containing all of the required inputs and outputs for the PV simulation
	std::unique_ptr<PVLossManager> m_PVLossManager;		/// An object containing the methods and data to construct the loss diagram
	//std::unique_ptr<PVSystemController> m_PVSystemController; /// An object that oversees the control of the PVSystem and any other technologies
};

/**
* \class PVLossManager
*
*  PVLossManager parses the PVIOManager and tracks detailed losses required for the loss diagram
*
*/
class PVLossManager
{
public:
	
	/// Construct a PVLossManager with a PVIOManager
	PVLossManager(PVIOManager * pvIOManager) :
		m_pvIOManager(pvIOManager){}

	/// Run the loss manager and parse the outputs
	bool Run();

private:

	/// A weak pointer to the PVIOManager, managed by PVSimulationManager
	PVIOManager * m_pvIOManager;
};
#endif
