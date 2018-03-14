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
	PVSimulationManager(std::shared_ptr<PVIOManager> IOManager, std::shared_ptr <PVLossManager> LossManager);

	/// PVSimulationManager is intended to be non-copyable
	PVSimulationManager(const PVSimulationManager&);
	PVSimulationManager& operator=(const PVSimulationManager&);

	/// Run a PV Simulation
	const bool Simulate();

private:

	/// Run a single step of the simulation
	const bool RunSingleStep();

	/// The member Manager objects are shared (though could be reconcieved as unique_ptr
	std::shared_ptr<PVIOManager> m_PVIOManager;
	std::shared_ptr<PVLossManager> m_PVLossManager;

	// SimulationManager uniquely manages ownership
	//std::unique_ptr<PVSystemController> m_PVSystemController;
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
	PVLossManager(std::shared_ptr<PVIOManager> pvIOManager) :
		m_pvIOManager(pvIOManager){}

	/// Run the loss manager and parse the outputs
	bool Run();

private:

	/// A shared pointer to the PVIOManager
	std::shared_ptr<PVIOManager> m_pvIOManager;
};
#endif
