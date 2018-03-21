#include "lib_pv_simulation_manager.h"

PVSimulationManager::PVSimulationManager(compute_module &cm)
{
	std::unique_ptr<PVIOManager> tmp(new PVIOManager(cm));
	m_PVIOManager = std::move(tmp);

	std::unique_ptr<PVLossManager> tmp2(new PVLossManager(m_PVIOManager.get()));
	m_PVLossManager = std::move(tmp2);

	std::unique_ptr<PVSystemController> tmp3(new PVSystemController(m_PVIOManager.get()));
	m_PVSystemController = std::move(tmp3);
}

const bool PVSimulationManager::Simulate()
{
	size_t index = 0;
	for (size_t year = 0; year != m_simulationIO->numberOfYears; year++) {
		for (size_t hour = 0; hour != util::hours_per_year; hour++) {
			for (size_t step = 0; step != m_simulationIO->stepsPerHour; step++) 
			{
				m_PVSystemController->RunSingleStep(index);
				index++;
			}

		}

	}

	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

const bool PVSimulationManager::RunSingleStep(size_t runIndex)
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}
