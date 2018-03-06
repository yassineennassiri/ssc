#ifndef _LIB_PV_SIMULATION_MANAGER_H_
#define _LIB_PV_SIMULATION_MANAGER_H_

#include <memory>

#include "lib_pv_io.h"

class PVLossManager;
class PVSystemController;

/** PVSimulationManager manages PV simulation details
* \class PVSimulationManager
*
* Given the PVIOManager and LossManager contructs and runs the SystemController
*/
class PVSimulationManager
{
public: 

	PVSimulationManager(std::shared_ptr<PVIOManager> IOManager, std::shared_ptr <PVLossManager> LossManager);

	// Can't copy a SimulationManager
	PVSimulationManager(const PVSimulationManager&);
	PVSimulationManager& operator=(const PVSimulationManager&);

	const bool Simulate();

private:

	const bool RunSingleStep();

	// Shared ownership
	std::shared_ptr<PVIOManager> m_PVIOManager;
	std::shared_ptr<PVLossManager> m_PVLossManager;

	// SimulationManager uniquely manages ownership
	//std::unique_ptr<PVSystemController> m_PVSystemController;
};

class PVLossManager
{
public:
	
	PVLossManager(std::shared_ptr<PVIOManager> pvIOManager) :
		m_pvIOManager(pvIOManager){}

	bool Run();

private:

	std::shared_ptr<PVIOManager> m_pvIOManager;
};
#endif
