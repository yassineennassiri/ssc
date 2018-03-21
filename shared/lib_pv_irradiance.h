#ifndef _LIB_PV_IRRADIANCE_H_
#define _LIB_PV_IRRADIANCE_H_

#include "lib_pv_io.h"
#include <memory>

// Forward declarations
class TrackerModel;
class ExternalShadeModel;
class SelfShadeModel;

/**
* \class IrradianceModel
*
* A IrradianceModel contains information relevent to one subarray, including which irradiance model should be used,
* and information about the irradiance properties for the subarray.
*
*/
class IrradianceModel
{
public:
	/// Construct a IrradianceModel from a PVIOManager
	IrradianceModel(PVIOManager * pvIOManager);

	/// Process the irradiance for one time step
	const bool RunSingleStep(const size_t runIndex);

private:

	/// Weak pointers that are owned by the PVSimulationManager or PVIOManager
	PVIOManager * m_pvIOManager;
	Irradiance_IO * m_irradianceIO;

	/// The IrradianceModel uniquely owns and manages these models
	std::unique_ptr<TrackerModel> m_TrackerModel;				/// The Tracker Model, which runs the tracking system for the subarray and modifies the POA irradiance
	std::unique_ptr<ExternalShadeModel> m_ExternalShadeModel;	/// The External Shade Model which accounts for shading of the subarray from objects, and modifies the POA irradiance and provides a DC loss component
	std::unique_ptr<SelfShadeModel> m_SelfShadeModel;			/// The Self Shade Model which accounts for shading of the subarray from other rows, and modifies the POA irradiance and provides a DC loss component
};

/**
* \class TrackerModel
*
* A TrackerModel contains the tracking information and model for one subarray
*
*/
class TrackerModel
{
public:
	/// Construct a TrackerModel from a PVIOManager
	TrackerModel(PVIOManager * pvIOManager);

	/// Process the irradiance for one time step
	const bool RunSingleStep(const size_t runIndex);

private:

	/// The PVIOManager is a weak pointer that is owned by the PVSimulationManager
	PVIOManager * m_pvIOManager;
};

/**
* \class ExternalShadeModel
*
* A ExternalShadeModel contains the near shading information and model for one subarray
*
*/
class ExternalShadeModel
{
public:
	/// Construct a ExternalShadeModel from a PVIOManager
	ExternalShadeModel(PVIOManager * pvIOManager);

	/// Process the irradiance for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a weak pointer that is owned by the PVSimulationManager
	PVIOManager * m_pvIOManager;
};

/**
* \class SelfShadeModel
*
* A SelfShadeModel contains the information about shading due to rows for one subarray
*
*/
class SelfShadeModel
{
public:
	/// Construct a SelfShadeModel from a PVIOManager
	SelfShadeModel(PVIOManager * pvIOManager);

	/// Process the irradiance for one time step
	const bool RunSingleStep();

private:

	/// The PVIOManager is a weak pointer that is owned by the PVSimulationManager
	PVIOManager * m_pvIOManager;
};

#endif