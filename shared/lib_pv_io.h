#ifndef _LIB_PV_IO_H_
#define _LIB_PV_IO_H_

#include <memory>
#include <math.h>

#include "lib_util.h"
#include "lib_irradproc.h"
#include "lib_weatherfile.h"

#include "../ssc/common.h"
#include "../ssc/core.h"

struct Simulation_IO;
struct Irradiance_IO;
struct Subarray_IO;

/*
struct PVSystem_IO;
struct MPPTController_IO;

struct Battery_IO;
struct Inverter_IO;
*/

class PVIOManager
{
public:

	PVIOManager(compute_module &cm);

	Simulation_IO * getSimulationIO() const;
	Irradiance_IO * getIrradianceIO() const;
	Subarray_IO * getSubarrayIO(size_t subarray) const;
	
	/*
	PVSystem_IO * getPVSystemIO() const;
	MPPTController_IO * getMPPTControllerIO() const;
	Battery_IO * getBatteryIO() const;
	Inverter_IO * getInverterIO() const;
	*/
	
private:
	
	// IOManager uniquely manages ownership
	std::unique_ptr<Irradiance_IO> m_SimulationIO;
	std::unique_ptr<Irradiance_IO> m_IrradianceIO;
	std::vector<std::unique_ptr<Subarray_IO>> m_SubarraysIO;

	/*
	std::unique_ptr<PVSystem_IO> m_PVSystemIO;
	std::unique_ptr<MPPTController_IO> m_MPPTControllerIO;
	std::unique_ptr<Battery_IO> m_BatteryIO;
	std::unique_ptr<Inverter_IO> m_InverterIO;
	*/

};


struct Subarray_IO
{
	Subarray_IO(compute_module &cm, size_t subarrayNumber)
	{
		std::string prefix = "subarray" + util::to_string(static_cast<int>(subarrayNumber));
		enable = cm.as_boolean(prefix + "enable");

		if (enable)
		{
			tilt = fabs(cm.as_double(prefix + "tilt"));
			azimuth = cm.as_double(prefix + "azimuth");
			trackMode = cm.as_integer(prefix + "track_mode");
			trackerRotationLimit = cm.as_double(prefix + "rotlim");
			tiltEqualLatitude = cm.as_boolean(prefix + "tile_eq_lat");
			groundCoverageRatio = cm.as_double(prefix + "gcr");
			monthlyTilt = cm.as_doublevec(prefix + "monthly_tilt");
			backtrackingEnabled = cm.as_boolean(prefix + "backtrack");
			derate = (1 - cm.as_double(prefix + "mismatch_loss") / 100) *
				(1 - cm.as_double(prefix + "diodeconn_loss") / 100) *
				(1 - cm.as_double(prefix + "dcwiring_loss") / 100) *
				(1 - cm.as_double(prefix + "tracking_loss") / 100) *
				(1 - cm.as_double(prefix + "nameplate_loss") / 100) *
				(1 - cm.as_double("dcoptimizer_loss") / 100);

			if (groundCoverageRatio < 0.01)
				throw compute_module::exec_error("pvsamv2", "array ground coverage ratio must obey 0.01 < gcr");
		}
	}

	/**
	* Subarray Inputs
	*/
	bool enable;
	size_t nStrings;
	std::vector<double> monthlySoiling;
	double derate;
	double dcLoss;
	double groundCoverageRatio;
	double tilt;
	double azimuth;
	int trackMode;
	double trackerRotationLimit;
	bool tiltEqualLatitude;
	std::vector<double> monthlyTilt;
	bool backtrackingEnabled;
	int shadeMode;

	/**
	*  Structure to represent the plane of array calculations returned by irradiance processor
	*/
	struct {
		double ibeam;
		double iskydiff;
		double ignddiff;
		double ipoa;
		int sunup;
		double aoi;
		double stilt;
		double sazi;
		double nonlinear_dc_shading_derate;
		bool usePOAFromWF;
		int poaShadWarningCount;
		poaDecompReq poaAll;
	} poa;

	/**
	* Structure to represent module level information, calculated by module model
	*/
	struct {
		double dcpwr;
		double dcv;
		double voc;
		double isc;
		double dceff;
		double tcell;
	} module;

};

struct Irradiance_IO 
{
	Irradiance_IO(compute_module &cm)
	{
		if (cm.is_assigned("solar_resource_file")){ 
			weatherDataProvider = std::unique_ptr<weather_data_provider>(new weatherfile(cm.as_string("solar_resource_file")));
		}
		else {
			weatherDataProvider = std::unique_ptr<weather_data_provider>(new weatherdata(cm.lookup("solar_resource_data")));
		}

		// Check weather file
		if (weatherDataProvider->has_message()) cm.log(weatherDataProvider->message(), SSC_WARNING);
		weatherfile *weatherFile = dynamic_cast<weatherfile*>(weatherDataProvider.get());
		if (!weatherFile->ok()) throw compute_module::exec_error("pvsamv2", weatherFile->message());
		if (weatherFile->has_message()) cm.log(weatherFile->message(), SSC_WARNING);
	}

	std::unique_ptr<weather_data_provider> weatherDataProvider;
};

struct Simulation_IO 
{
	Simulation_IO(compute_module &cm, Irradiance_IO & IrradianceIO)
	{
		numberOfWeatherFileRecords = IrradianceIO.weatherDataProvider->nrecords();
		stepsPerHour = numberOfWeatherFileRecords / 8760;
		dtHour = 1.0 / stepsPerHour;
		useLifetimeOutput = cm.as_integer("system_use_lifetime_output");
		numberOfYears = 1;
		if (useLifetimeOutput){
			numberOfYears = cm.as_integer("analysis_period");
		}
	}

	size_t numberOfYears;
	size_t numberOfWeatherFileRecords;
	size_t stepsPerHour;
	double dtHour;
	bool useLifetimeOutput;
};

#endif