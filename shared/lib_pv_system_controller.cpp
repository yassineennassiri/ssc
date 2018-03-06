#include "lib_pv_system_controller.h"

PVSystemController::PVSystemController(std::shared_ptr <PVIOManager> pvIOManager) :
	m_pvIOManager(pvIOManager)
{
	std::unique_ptr<PVSystem> ptr(new PVSystem(pvIOManager));
	m_pvSystem = std::move(ptr);
}

const bool PVSystemController::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

PVSystem::PVSystem(std::shared_ptr<PVIOManager> pvIOManager) :
	m_pvIOManager(pvIOManager)
{

}

const bool PVSystem::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

