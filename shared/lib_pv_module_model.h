#ifndef _LIB_PV_MODULE_MODEL_
#define _LIB_PV_MODULE_MODEL_

#include <memory>

// Forward declarations
class SpectralModel;
class IAMModel;
class ModuleThermalModel;
class ModuleCoverModel;

/**
* \class ModuleModel
*
* A ModuleModel contains information relevent to one subarray, including which type of module is installed, and
* what models should be used to evaluate performance
*
*/
class ModuleModel
{
public:
	/// Construct a ModuleModel from a PVIOManager
	ModuleModel(PVIOManager * pvIOManager);

	/// Process the irradiance for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a weak pointer that is owned by the PVSimulationManager
	PVIOManager * m_pvIOManager;

	/// The ModuleModel uniquely owns and manages these models
	std::unique_ptr<SpectralModel> m_spectralModel;				/// The Spectral model, which modifies the DC power depending on the spectrum of light energy
	std::unique_ptr<IAMModel> m_IAMModel;						/// The IAM (incident angle modifier) model, which modifies the DC power depending on the angle of incidence
	std::unique_ptr<ModuleThermalModel> m_moduleThermalModel;	/// The Module Thermal model, which modifies the DC power depending on the module temperature 
	std::unique_ptr<ModuleCoverModel> m_moduleCoverModel;		/// The Cover Model which modifies the DC power depending on the module cover
};

/**
* \class IAMModel
*
* A IAMModel contains the incidence angle modifier information and model for one subarray
*
*/
class IAMModel
{
public:
	/// Construct a IAMModel from a PVIOManager
	IAMModel(PVIOManager * pvIOManager);

	/// Process the irradiance for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a weak pointer that is owned by the PVSimulationManager
	PVIOManager * m_pvIOManager;
};

/**
* \class SpectralModel
*
* A SpectralModel contains the spectral information and model for one subarray
*
*/
class SpectralModel
{
public:
	/// Construct a ModuleModel from a PVIOManager
	SpectralModel(PVIOManager * pvIOManager);

	/// Process the irradiance for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a weak pointer that is owned by the PVSimulationManager
	PVIOManager * m_pvIOManager;
};

/**
* \class ModuleThermalModel
*
* A ModuleThermalModel contains the thermal properties and model for one subarray
*
*/
class ModuleThermalModel
{
public:
	/// Construct a ModuleThermalModel from a PVIOManager
	ModuleThermalModel(PVIOManager * pvIOManager);

	/// Process the irradiance for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a weak pointer that is owned by the PVSimulationManager
	PVIOManager * m_pvIOManager;
};

/**
* \class ModuleCoverModel
*
* A ModuleCoverModel contains the cover information and model for one subarray
*
*/
class ModuleCoverModel
{
public:
	/// Construct a ModuleThermalModel from a PVIOManager
	ModuleCoverModel(PVIOManager * pvIOManager);

	/// Process the irradiance for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a weak pointer that is owned by the PVSimulationManager
	PVIOManager * m_pvIOManager;
};

#endif