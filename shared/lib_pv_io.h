#ifndef _LIB_PV_IO_H_
#define _LIB_PV_IO_H_

#include <map>
#include <memory>
#include <math.h>

#include "lib_util.h"
#include "lib_irradproc.h"
#include "lib_weatherfile.h"

#include "../ssc/common.h"
#include "../ssc/core.h"

/// Structure containing data relevent at the SimulationManager level
struct Simulation_IO;

/// Structure contain data relevent to the Irradiance model
struct Irradiance_IO;

/// Structure containing subarray-level IO information
struct Subarray_IO;

/// Structure containing inverter IO information
struct Inverter_IO;

/// Structure containing MPPT controller IO information
struct MPPTController_IO;

/*
struct PVSystem_IO;
struct Battery_IO;
*/

/**
* \class PVIOManager
*
* This class contains the input and output data needed by all of the submodels in the detailed PV model
* It is intended to be passed as a pointer and modified with the goal of encapsulating all input that comes
* in at the user level and all output that needs to ultimately be passed back to the user.  
*
*
* \note The PVIOManager contains structures that contain specific information about each system component
*/
class PVIOManager
{
public:
	
	/// Create a PVIOManager object by parsing the compute model
	PVIOManager(compute_module &cm);

	/// Return Simulation specific information
	Simulation_IO * getSimulationIO() const;

	/// Return Irradiance specific information
	Irradiance_IO * getIrradianceIO() const;

	/// Return Subarray specific information for the given subarray
	Subarray_IO * getSubarrayIO(size_t subarray) const;

	/// Return MPPT Controller specific information for the given subarray
	MPPTController_IO * getMPPTControllerIO() const;
	
	
	/*
	PVSystem_IO * getPVSystemIO() const;
	Battery_IO * getBatteryIO() const;
	Inverter_IO * getInverterIO() const;
	*/
	
private:
	
	/** These structures contain specific IO data for each part of the model
	  * They are owned exclusively by the PVIOManager 
	  */
	std::unique_ptr<Irradiance_IO> m_SimulationIO;
	std::unique_ptr<Irradiance_IO> m_IrradianceIO;
	std::vector<std::unique_ptr<Subarray_IO>> m_SubarraysIO;
	std::unique_ptr<MPPTController_IO> m_MPPTControllerIO;
	std::unique_ptr<Inverter_IO> m_InverterIO;

	/*
	std::unique_ptr<PVSystem_IO> m_PVSystemIO;
	std::unique_ptr<Battery_IO> m_BatteryIO;
	*/

};


struct Subarray_IO
{
	Subarray_IO(compute_module &cm, size_t subarrayNumber)
	{
		std::string prefix = "subarray" + util::to_string(static_cast<int>(subarrayNumber)) + "_";
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
		if (cm.is_assigned("solar_resource_file")) {
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
		if (useLifetimeOutput) {
			numberOfYears = cm.as_integer("analysis_period");
		}
	}

	size_t numberOfYears;
	size_t numberOfWeatherFileRecords;
	size_t stepsPerHour;
	double dtHour;
	bool useLifetimeOutput;
};

struct Inverter_IO
{
	Inverter_IO(compute_module &cm)
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

	int inverterType;
	double efficiency;
	double Paco;
	double Pdco;
	double Vdco;
	double Pso;
	double Pntare;
	double C0;
	double C1;
	double C2;
	double C3;
	double ratedACOuput;
	std::vector<double> partloadPowerPercent;
	std::vector<double> partloadEfficiency;
};

struct MPPTController_IO
{
	MPPTController_IO(compute_module &cm)
	{
		n_enabledSubarrays = 0;
		for (int subarrayNumber = 0; subarrayNumber != 4; subarrayNumber++)
		{
			std::string prefix = "subarray" + util::to_string(subarrayNumber) + "_";
			bool enable = cm.as_boolean(prefix + "enable");
			if (enable)
			{
				int MPPTType = cm.as_integer(prefix + "mppt_type");
				int MPPTPort = cm.as_integer(prefix + "mppt_port");
				subarrayMPPTControllers[subarrayNumber] = MPPTType;
				subarrayMPPTPorts[subarrayNumber] = MPPTPort;
				n_enabledSubarrays++;
			}
		}
	}
	int n_enabledSubarrays;
	std::map<const int, int > subarrayMPPTControllers;
	std::map<const int, int > subarrayMPPTPorts;
};


#endif