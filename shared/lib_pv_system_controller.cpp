#include "lib_pv_system_controller.h"

PVSystemController::PVSystemController(PVIOManager * pvIOManager) :
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

PVSystem::PVSystem(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
	std::unique_ptr<PVDCController> ptr(new PVDCController(pvIOManager));
	m_pvDCController = std::move(ptr);

	std::unique_ptr<PVACController> ptr2(new PVACController(pvIOManager));
	m_pvACController = std::move(ptr2);
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

PVDCController::PVDCController(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

PVACController::PVACController(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

MPPTController::MPPTController(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
	m_MPPTIO = pvIOManager->getMPPTControllerIO();
	for (size_t subarrayNumber = 0; subarrayNumber != m_MPPTIO->n_enabledSubarrays; subarrayNumber++)
	{
		std::unique_ptr<Subarray> ptr(new Subarray(pvIOManager));
		m_Subarrays.push_back(std::move(ptr));
	}
}

const bool MPPTController::RunSingleStep()
{
	for (size_t subarrayNumber = 0; subarrayNumber != m_MPPTIO->n_enabledSubarrays; subarrayNumber++)
	{
		m_Subarrays[subarrayNumber]->RunSingleStep();
		
		// do stuff with MPPT feedback to subarray and check error conditions
	}
	return EXIT_SUCCESS;
}

Subarray::Subarray(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool Subarray::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}
