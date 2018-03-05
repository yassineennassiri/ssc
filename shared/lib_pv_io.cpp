#include <memory>
#include <vector>

#include "lib_pv_io.h"

PVIOManager::PVIOManager(compute_module & cm)
{
	for (size_t subarray = 0; subarray != 4; subarray++)
	{
		std::unique_ptr<Subarray_IO> ptr(new Subarray_IO(cm, subarray));
		m_SubarraysIO.push_back(std::move(ptr));
	}
}