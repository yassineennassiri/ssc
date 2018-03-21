#include "lib_pv_irradiance.h"

IrradianceModel::IrradianceModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager), 
	m_irradianceIO(pvIOManager->getIrradianceIO())
{
}

const bool IrradianceModel::RunSingleStep(const size_t runIndex)
{

	if (!m_irradianceIO->weatherDataProvider.get()->read(&m_irradianceIO->weatherRecord)) {
		throw compute_module::exec_error("pvsamv2", "could not read data line " + util::to_string((int)(runIndex + 1)) + " in weather file while loading POA data");
	}

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

const bool TrackerModel::RunSingleStep(const size_t runIndex)
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