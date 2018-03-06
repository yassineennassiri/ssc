#ifndef _LIB_PV_SYSTEM_CONTROLLER_H_
#define _LIB_PV_SYSTEM_CONTROLLER_H_

#include <memory>

#include "lib_pv_io.h"

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
};

#endif
