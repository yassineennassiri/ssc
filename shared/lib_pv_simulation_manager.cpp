#include "lib_pv_simulation_manager.h"

PVSimulationManager::PVSimulationManager(compute_module &cm)
{
	std::unique_ptr<PVIOManager> tmp(new PVIOManager(cm));
	m_PVIOManager = std::move(tmp);

	std::unique_ptr<PVLossManager> tmp2(new PVLossManager(m_PVIOManager.get()));
	m_PVLossManager = std::move(tmp2);
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
