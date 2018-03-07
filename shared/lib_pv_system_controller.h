#ifndef _LIB_PV_SYSTEM_CONTROLLER_H_
#define _LIB_PV_SYSTEM_CONTROLLER_H_

#include <memory>

#include "lib_pv_io.h"

class Subarray;
class MPPTController;
class PVDCController;
class PVACController;
class PVSystem;
class PVSystemController
{
public:
	PVSystemController(std::shared_ptr<PVIOManager> pvIOManager);

	const bool RunSingleStep();

private:
	std::shared_ptr<PVIOManager> m_pvIOManager;
	std::unique_ptr<PVSystem> m_pvSystem;
};

class PVSystem
{
public:
	PVSystem(std::shared_ptr<PVIOManager> pvIOManager);

	const bool RunSingleStep();

private:
	std::shared_ptr<PVIOManager> m_pvIOManager;
	std::unique_ptr<PVDCController> m_pvDCController;
	std::unique_ptr<PVACController> m_pvACController;
};

class PVDCController
{
public: 
	PVDCController(std::shared_ptr<PVIOManager> pvIOManager);

	//const bool RunSingleStep();
private:
	std::shared_ptr<PVIOManager> m_pvIOManager;
	std::vector<std::unique_ptr<MPPTController>> m_MPPTControllers;

};

class PVACController
{
public:
	PVACController(std::shared_ptr<PVIOManager> pvIOManager);

	//const bool RunSingleStep();

private:
	std::shared_ptr<PVIOManager> m_pvIOManager;

};

class MPPTController
{
public:
	MPPTController(std::shared_ptr<PVIOManager> pvIOManager);

	const bool RunSingleStep();

private:
	std::shared_ptr<PVIOManager> m_pvIOManager;
	std::vector<std::unique_ptr<Subarray>> m_Subarrays;

	MPPTController_IO * m_MPPTIO;
};

class Subarray
{
public:
	Subarray(std::shared_ptr<PVIOManager> pvIOManager);

	const bool RunSingleStep();

private:
	std::shared_ptr<PVIOManager> m_pvIOManager;
	
};

#endif
