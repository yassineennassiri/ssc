#include "lib_pv_irradiance.h"

IrradianceModel::IrradianceModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool IrradianceModel::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

TrackerModel::TrackerModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool TrackerModel::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

ExternalShadeModel::ExternalShadeModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool ExternalShadeModel::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}

SelfShadeModel::SelfShadeModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager)
{
}

const bool SelfShadeModel::RunSingleStep()
{
	// do stuff
	bool ranSuccessfully = true;

	if (ranSuccessfully)
		return EXIT_SUCCESS;
	else
		return EXIT_FAILURE;
}