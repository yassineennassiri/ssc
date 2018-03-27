#include "lib_pv_irradiance.h"

IrradianceModel::IrradianceModel(PVIOManager * pvIOManager) :
	m_pvIOManager(pvIOManager), 
	m_irradianceIO(pvIOManager->getIrradianceIO())
{
	// register the subarrays with the IrradianceModel
	for (size_t subarrayNumber = 0; subarrayNumber != m_pvIOManager->getMPPTControllerIO()->maxSubarrays; subarrayNumber++){
		m_subarraysIO.push_back(m_pvIOManager->getSubarrayIO(subarrayNumber));
	}
}

const bool IrradianceModel::RunSingleStep(const size_t runIndex)
{
	// convenience variables 
	compute_module * computeModule = m_pvIOManager->getComputeModule();
	weather_record wf = m_irradianceIO->weatherRecord;
	weather_header hdr = m_irradianceIO->weatherHeader;
	int radiationMode = m_irradianceIO->radiationMode;
	int skyModel = m_irradianceIO->skyModel;
	int irradiationMax = m_irradianceIO->irradiationMax;
	bool instantaneous = m_irradianceIO->instantaneous;
	double dtHour = m_irradianceIO->dtHour;
	bool useWeatherFileAlbedo = m_irradianceIO->useWeatherFileAlbedo;
	std::vector<double> monthlyAlbedo = m_irradianceIO->userSpecifiedMonthlyAlbedo;
	double albedo = 0.2;

	// need to define template functions for util::to_string, and compute_module::log
	float runIndexF = static_cast<float>(runIndex);
	int runIndexI = static_cast<int>(runIndex);

	// Read one line from the weather file
	if (!m_irradianceIO->weatherDataProvider.get()->read(&wf)) {
		throw compute_module::exec_error("pvsamv2", "could not read data line " + util::to_string(runIndexI+1) + " in weather file");
	}

	if (useWeatherFileAlbedo && std::isfinite(wf.alb) && wf.alb > 0 && wf.alb < 1)
		albedo = wf.alb;
	else if (wf.month >= 1 && wf.month <= 12)
		albedo = monthlyAlbedo[wf.month - 1];
	else
		throw compute_module::exec_error("pvsamv2",
			util::format("Error retrieving albedo value: Invalid month in weather file or invalid albedo value in weather file"));

	// Check for missing data
	if ((wf.gh != wf.gh) && (radiationMode == DN_GH || radiationMode == GH_DF)) {
		computeModule->log(util::format("missing global irradiance %lg W/m2 at time [y:%d m:%d d:%d h:%d], exiting",
			wf.gh, wf.year, wf.month, wf.day, wf.hour), SSC_ERROR, runIndexF);
		return EXIT_FAILURE;
	}
	if ((wf.dn != wf.dn) && (radiationMode == DN_DF || radiationMode == DN_GH)) {
		computeModule->log(util::format("missing beam irradiance %lg W/m2 at time [y:%d m:%d d:%d h:%d], exiting",
			wf.dn, wf.year, wf.month, wf.day, wf.hour), SSC_ERROR, runIndexF);
		return EXIT_FAILURE;
	}
	if ((wf.df != wf.df) && (radiationMode == DN_DF || radiationMode == GH_DF)) {
		computeModule->log(util::format("missing diffuse irradiance %lg W/m2 at time [y:%d m:%d d:%d h:%d], exiting",
			wf.df, wf.year, wf.month, wf.day, wf.hour), SSC_ERROR, runIndexF);
		return EXIT_FAILURE;
	}
	if ((wf.poa != wf.poa) && (radiationMode == POA_R || radiationMode == POA_P)) {
		computeModule->log(util::format("missing POA irradiance %lg W/m2 at time [y:%d m:%d d:%d h:%d], exiting",
			wf.poa, wf.year, wf.month, wf.day, wf.hour), SSC_ERROR, runIndexF);
		return EXIT_FAILURE;
	}
	if (wf.tdry != wf.tdry) {
		computeModule->log(util::format("missing temperature %lg W/m2 at time [y:%d m:%d d:%d h:%d], exiting",
			wf.tdry, wf.year, wf.month, wf.day, wf.hour), SSC_ERROR, runIndexF);
		return EXIT_FAILURE;
	}
	if (wf.wspd != wf.wspd) {
		computeModule->log(util::format("missing wind speed %lg W/m2 at time [y:%d m:%d d:%d h:%d], exiting",
			wf.wspd, wf.year, wf.month, wf.day, wf.hour), SSC_ERROR, runIndexF);
		return EXIT_FAILURE;
	}

	// Check for bad data
	if ((wf.gh < 0 || wf.gh > irradiationMax) && (radiationMode == DN_GH || radiationMode == GH_DF))
	{
		computeModule->log(util::format("out of range global irradiance %lg W/m2 at time [y:%d m:%d d:%d h:%d], set to zero",
			wf.gh, wf.year, wf.month, wf.day, wf.hour), SSC_WARNING, runIndexF);
		wf.gh = 0;
	}
	if ((wf.dn < 0 || wf.dn > irradiationMax) && (radiationMode == DN_DF || radiationMode == DN_GH))
	{
		computeModule->log(util::format("out of range beam irradiance %lg W/m2 at time [y:%d m:%d d:%d h:%d], set to zero",
			wf.dn, wf.year, wf.month, wf.day, wf.hour), SSC_WARNING, runIndexF);
		wf.dn = 0;
	}
	if ((wf.df < 0 || wf.df > irradiationMax) && (radiationMode == DN_DF || radiationMode == GH_DF))
	{
		computeModule->log(util::format("out of range diffuse irradiance %lg W/m2 at time [y:%d m:%d d:%d h:%d], set to zero",
			wf.df, wf.year, wf.month, wf.day, wf.hour), SSC_WARNING, runIndexF);
		wf.df = 0;
	}
	if ((wf.poa < 0 || wf.poa > irradiationMax) && (radiationMode == POA_R || radiationMode == POA_P))
	{
		computeModule->log(util::format("out of range POA irradiance %lg W/m2 at time [y:%d m:%d d:%d h:%d], set to zero",
			wf.poa, wf.year, wf.month, wf.day, wf.hour), SSC_WARNING, runIndexF);
		wf.poa = 0;
	}

	/*  For this piece will require:
	1. Having calculated the components of radiation for input POA
	2. Interacting with the Subarray_IO for the tracking and orientation info required to calculate POA irradiance on each subarray

	irrad irr;
	irr.set_time(wf.year, wf.month, wf.day, wf.hour, wf.minute,
		instantaneous ? IRRADPROC_NO_INTERPOLATE_SUNRISE_SUNSET : dtHour);
	irr.set_location(hdr.lat, hdr.lon, hdr.tz);

	irr.set_sky_model(skyModel, albedo);
	if (radiationMode == DN_DF) irr.set_beam_diffuse(wf.dn, wf.df);
	else if (radiationMode == DN_GH) irr.set_global_beam(wf.gh, wf.dn);
	else if (radiationMode == GH_DF) irr.set_global_diffuse(wf.gh, wf.df);

	

	int returnCode = irr.calc();
	if (returnCode != 0)
		throw compute_module::exec_error("pvsamv2",
			util::format("failed to calculate irradiance incident on surface (POA) %d (code: %d) [y:%d m:%d d:%d h:%d]",
				nn + 1, code, wf.year, wf.month, wf.day, wf.hour));

	*/

	// if year one.  There are still other outputs
	m_irradianceIO->p_weatherFileGHI[runIndex] = static_cast<ssc_number_t>(wf.gh);
	m_irradianceIO->p_weatherFileDNI[runIndex] = static_cast<ssc_number_t>(wf.dn);
	m_irradianceIO->p_weatherFileDHI[runIndex] = static_cast<ssc_number_t>(wf.df);
	m_irradianceIO->p_weatherFileWindSpeed[runIndex] = static_cast<ssc_number_t>(wf.wspd);
	m_irradianceIO->p_weatherFileAmbientTemp[runIndex] = static_cast<ssc_number_t>(wf.tdry);
	m_irradianceIO->p_weatherFileAlbedo[runIndex] = static_cast<ssc_number_t>(albedo);
	m_irradianceIO->p_weatherFileSnowDepth[runIndex] = static_cast<ssc_number_t>(wf.snow);

	// if made it here, return success
	return EXIT_SUCCESS;
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