#include "lib_pv_io.h"
#include "lib_pv_module_model.h"

ModuleModel::ModuleModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool ModuleModel::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

SpectralModel::SpectralModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool SpectralModel::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

IAMModel::IAMModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool IAMModel::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

ModuleThermalModel::ModuleThermalModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool ModuleThermalModel::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

ModuleCoverModel::ModuleCoverModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool ModuleCoverModel::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}