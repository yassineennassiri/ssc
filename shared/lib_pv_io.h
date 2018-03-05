#ifndef _LIB_PV_IO_H_
#define _LIB_PV_IO_H_

#include <memory>
#include "../ssc/core.h"


struct PVSystem_IO;
struct MPPTController_IO;
struct Subarray_IO;
struct Irradiance_IO;
struct Module_IO;
struct Tracker_IO;
struct Shade_IO;
struct Spectral_IO;
struct IAM_IO;
struct Battery_IO;
struct Inverter_IO;
struct InverterThermal_IO;

class PVIOManager
{
public:

	PVIOManager(compute_module &cm);

	PVIOManager(const PVIOManager&);
	PVIOManager& operator=(const PVIOManager&);

	PVSystem_IO * getPVSystemIO() const;
	MPPTController_IO * getMPPTControllerIO() const;
	Subarray_IO * getSubarrayIO() const;
	Irradiance_IO * getIrradianceIO() const;
	Module_IO * getModuleIO() const;
	Tracker_IO * getTrackerIO() const;
	Shade_IO * getShadeIO() const;
	Spectral_IO * getSpectralIO() const;
	IAM_IO * getIAMIO() const;
	Battery_IO * getBatteryIO() const;
	Inverter_IO * getInverterIO() const;
	InverterThermal_IO * getInverterThermalIO() const;


private:

	// IOManager uniquely manages ownership
	std::unique_ptr<PVSystem_IO> m_PVSystemIO;
	std::unique_ptr<MPPTController_IO> m_MPPTControllerIO;
	std::unique_ptr<Subarray_IO> m_SubarrayIO;
	std::unique_ptr<Irradiance_IO> m_IrradianceIO;
	std::unique_ptr<Module_IO> m_ModuleIO;
	std::unique_ptr<Tracker_IO> m_TrackerIO;
	std::unique_ptr<Shade_IO> m_ShadeIO;
	std::unique_ptr<Spectral_IO> m_SpectralIO;
	std::unique_ptr<IAM_IO> m_ModuleIO;
	std::unique_ptr<Battery_IO> m_BatteryIO;
	std::unique_ptr<Inverter_IO> m_InverterIO;
	std::unique_ptr<InverterThermal_IO> m_InverterThermal_IO;

};

#endif