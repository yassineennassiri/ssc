#include <memory>
#include <vector>

#include "lib_pv_io.h"

PVIOManager::PVIOManager(compute_module & cm)
{
	std::unique_ptr<Irradiance_IO> ptr(new Irradiance_IO(cm));
	m_IrradianceIO = std::move(ptr);

	std::unique_ptr<Simulation_IO> ptr2(new Simulation_IO(cm, *m_IrradianceIO));
	m_SimulationIO = std::move(ptr2);

	std::unique_ptr<MPPTController_IO> ptr3(new MPPTController_IO(cm));
	m_MPPTControllerIO = std::move(ptr3);

	for (size_t subarray = 1; subarray <= m_MPPTControllerIO->maxSubarrays; subarray++)
	{
		std::unique_ptr<Subarray_IO> ptr4(new Subarray_IO(cm, subarray));
		m_SubarraysIO.push_back(std::move(ptr4));
	}
	m_computeModule = &cm;
}
Irradiance_IO * PVIOManager::getIrradianceIO() const { return m_IrradianceIO.get(); }
MPPTController_IO * PVIOManager::getMPPTControllerIO() const {return m_MPPTControllerIO.get();}
compute_module * PVIOManager::getComputeModule() const { return m_computeModule; }
Subarray_IO * PVIOManager::getSubarrayIO(size_t subarrayNumber) const { return m_SubarraysIO[subarrayNumber].get(); }
Simulation_IO * PVIOManager::getSimulationIO() const { return m_SimulationIO.get(); }


Irradiance_IO::Irradiance_IO(compute_module &cm)
{
	radiationMode = cm.as_integer("irrad_mode");
	skyModel = cm.as_integer("sky_model");

	if (cm.is_assigned("solar_resource_file")) {
		weatherDataProvider = std::unique_ptr<weather_data_provider>(new weatherfile(cm.as_string("solar_resource_file")));
	}
	else if (cm.is_assigned("solar_resource_data")) {
		weatherDataProvider = std::unique_ptr<weather_data_provider>(new weatherdata(cm.lookup("solar_resource_data")));
	}
	else {
		throw compute_module::exec_error("pvsamv2", "No weather data supplied");
	}

	// Check weather file
	if (weatherDataProvider->has_message()) cm.log(weatherDataProvider->message(), SSC_WARNING);
	weatherfile *weatherFile = dynamic_cast<weatherfile*>(weatherDataProvider.get());
	if (!weatherFile->ok()) throw compute_module::exec_error("pvsamv2", weatherFile->message());
	if (weatherFile->has_message()) cm.log(weatherFile->message(), SSC_WARNING);

	// assumes instantaneous values, unless hourly file with no minute column specified
	tsShiftHours = 0.0;
	instantaneous = true;
	if (weatherDataProvider->has_data_column(weather_data_provider::MINUTE))
	{
		// if we have an file with a minute column, then
		// the starting time offset equals the time 
		// of the first record (for correct plotting)
		// this holds true even for hourly data with a minute column
		weather_record rec;
		if (weatherDataProvider->read(&rec))
			tsShiftHours = rec.minute / 60.0;

		weatherDataProvider->rewind();
	}
	else if (weatherDataProvider->nrecords() == 8760)
	{
		// hourly file with no minute data column.  assume
		// integrated/averaged values and use mid point convention for interpreting results
		instantaneous = false;
		tsShiftHours = 0.5;
	}
	else
		throw compute_module::exec_error("pvsamv2", "subhourly weather files must specify the minute for each record");

	weatherDataProvider->header(&weatherHeader);

	//total number of records in the weather file (i.e. 8760 * timestep)
	numberOfWeatherFileRecords = weatherDataProvider->nrecords();
	stepsPerHour = numberOfWeatherFileRecords / 8760;
	dtHour = 1.0 / stepsPerHour;

	if (stepsPerHour < 1 || stepsPerHour > 60 || stepsPerHour * 8760 != numberOfWeatherFileRecords)
		throw compute_module::exec_error("pvsamv2", util::format("invalid number of data records (%zu): must be an integer multiple of 8760", numberOfWeatherFileRecords));

	useWeatherFileAlbedo = cm.as_boolean("use_wf_albedo");
	userSpecifiedMonthlyAlbedo = cm.as_doublevec("albedo");

	// Ensure all outputs are allocated
	AllocateOutputs(cm);
}

void Irradiance_IO::AllocateOutputs(compute_module &cm)
{
	p_weatherFileGHI = cm.allocate("gh", numberOfWeatherFileRecords);
	p_weatherFileDNI = cm.allocate("dn", numberOfWeatherFileRecords);
	p_weatherFileDHI = cm.allocate("df", numberOfWeatherFileRecords);
	p_sunPositionTime = cm.allocate("sunpos_hour", numberOfWeatherFileRecords);
	p_weatherFileWindSpeed = cm.allocate("wspd", numberOfWeatherFileRecords);
	p_weatherFileAmbientTemp = cm.allocate("tdry", numberOfWeatherFileRecords);
	p_weatherFileAlbedo = cm.allocate("alb", numberOfWeatherFileRecords);
	p_weatherFileSnowDepth = cm.allocate("snowdepth", numberOfWeatherFileRecords);

	// If using input POA, must have POA for every subarray or assume POA applies to each subarray
	if (radiationMode == POA_R || radiationMode == POA_P) {
		for (size_t subarray = 0; subarray != numberOfSubarrays; subarray++) {
			std::string wfpoa = "wfpoa" + util::to_string(static_cast<int>(subarray + 1));
			p_weatherFilePOA.push_back(cm.allocate(wfpoa, numberOfWeatherFileRecords));
		}
	}

	//set up the calculated components of irradiance such that they aren't reported if they aren't assigned
	//three possible calculated irradiance: gh, df, dn
	if (radiationMode == DN_DF) p_IrradianceCalculated[0] = cm.allocate("gh_calc", numberOfWeatherFileRecords); //don't calculate global for POA models
	if (radiationMode == DN_GH || radiationMode == POA_R || radiationMode == POA_P) p_IrradianceCalculated[1] = cm.allocate("df_calc", numberOfWeatherFileRecords);
	if (radiationMode == GH_DF || radiationMode == POA_R || radiationMode == POA_P) p_IrradianceCalculated[2] = cm.allocate("dn_calc", numberOfWeatherFileRecords);

	//output arrays for solar position calculations- same for all four subarrays
	p_sunZenithAngle = cm.allocate("sol_zen", numberOfWeatherFileRecords);
	p_sunAltitudeAngle = cm.allocate("sol_alt", numberOfWeatherFileRecords);
	p_sunAzimuthAngle = cm.allocate("sol_azi", numberOfWeatherFileRecords);
	p_absoluteAirmass = cm.allocate("airmass", numberOfWeatherFileRecords);
	p_sunUpOverHorizon = cm.allocate("sunup", numberOfWeatherFileRecords);
}

void Irradiance_IO::AssignOutputs(compute_module &cm)
{
	cm.assign("ts_shift_hours", var_data((ssc_number_t)tsShiftHours));
}

Simulation_IO::Simulation_IO(compute_module &cm, Irradiance_IO & IrradianceIO)
{
	numberOfWeatherFileRecords = IrradianceIO.numberOfWeatherFileRecords;
	stepsPerHour = IrradianceIO.stepsPerHour;
	dtHour = IrradianceIO.dtHour;

	useLifetimeOutput = cm.as_integer("system_use_lifetime_output");
	numberOfYears = 1;
	if (useLifetimeOutput) {
		numberOfYears = cm.as_integer("analysis_period");
	}
	numberOfSteps = numberOfYears * numberOfWeatherFileRecords;
}

Inverter_IO::Inverter_IO(compute_module &cm)
{
	inverterType = cm.as_integer("inverter_model");
	std::string type;
	if (inverterType == 0)
		type = "snl";
	else if (inverterType == 1)
		type = "ds";
	else if (inverterType == 2) {
		type = "pd";
		partloadPowerPercent = cm.as_doublevec("inv_pd_partload");
		partloadEfficiency = cm.as_doublevec("inv_pd_efficiency");
	}
	else if (inverterType == 3)
		type = "cec_cg";
	else
		throw compute_module::exec_error("pvsamv2", "invalid inverter model type");

	Paco = cm.as_double("inv_" + type + "_paco");
	Pdco = cm.as_double("inv_" + type + "_pdco");
	Vdco = cm.as_double("inv_" + type + "_vdco");
	Pso = cm.as_double("inv_" + type + "_pso");
	Pntare = cm.if_assigned_as_double("inv_" + type + "_pnt");
	C0 = C1 = C2 = 0;
	if (cm.is_assigned("inv_" + type + "_c0"))
		C0 = cm.as_double("inv_" + type + "c0");
	if (cm.is_assigned("inv_" + type + "_c1"))
		C1 = cm.as_double("inv_" + type + "c1");
	if (cm.is_assigned("inv_" + type + "_c2"))
		C2 = cm.as_double("inv_" + type + "c2");
	efficiency = cm.if_assigned_as_double("inv_" + type + "_eff");
	ratedACOuput = Paco;
}

MPPTController_IO::MPPTController_IO(compute_module &cm)
{
	n_enabledSubarrays = 0;
	size_t tmp_max = 0;
	for (size_t subarrayNumber = 1; subarrayNumber <= maxSubarrays; subarrayNumber++)
	{
		std::string prefix = "subarray" + util::to_string(static_cast<int>(subarrayNumber)) + "_";
		bool enable = true;
		if (subarrayNumber > 1) {
			enable = cm.as_boolean(prefix + "enable");
		}
		if (enable)
		{
			// need to define these in var table
			/*
			size_t MPPTType = static_cast<size_t>(cm.as_integer(prefix + "mppt_type"));
			size_t MPPTPort = static_cast<size_t>(cm.as_integer(prefix + "mppt_port"));
			*/
			size_t MPPTType = 1;
			size_t MPPTPort = 1;
			subarrayMPPTControllers[subarrayNumber] = MPPTType;
			subarrayMPPTPorts[subarrayNumber] = MPPTPort;
			n_enabledSubarrays++;

			if (MPPTType > tmp_max)
				tmp_max = MPPTType;
		}
	}
	n_MPPTControllers = tmp_max;
}