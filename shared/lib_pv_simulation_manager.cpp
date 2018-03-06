#include "lib_pv_simulation_manager.h"

PVSimulationManager::PVSimulationManager(std::shared_ptr<PVIOManager> IOManager, std::shared_ptr<PVLossManager> LossManager) :
	m_PVIOManager(IOManager),
	m_PVLossManager(LossManager)
{
	
}

const bool PVSimulationManager::Simulate()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

const bool PVSimulationManager::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}
