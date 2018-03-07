#include <memory>
#include <vector>

#include "lib_pv_io.h"

PVIOManager::PVIOManager(compute_module & cm)
{
	std::unique_ptr<Irradiance_IO> ptr(new Irradiance_IO(cm));
	m_IrradianceIO = std::move(ptr);

	for (size_t subarray = 0; subarray != 4; subarray++)
	{
		std::unique_ptr<Subarray_IO> ptr(new Subarray_IO(cm, subarray));
		m_SubarraysIO.push_back(std::move(ptr));
	}
}

MPPTController_IO * PVIOManager::getMPPTControllerIO() const
{
	return m_MPPTControllerIO.get();
}