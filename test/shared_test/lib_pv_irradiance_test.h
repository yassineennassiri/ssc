#ifndef _LIB_PV_IRRADIANCE_TEST_H_
#define _LIB_PV_IRRADIANCE_TEST_H_


#include <lib_pv_irradiance.h>
#include "../input_cases/pvsamv1_cases.h"
#include <cmod_pvsamv2.h>

class LibPVIrradiance : public ::testing::Test {

public:

	// Member variables, including PVIOManager and IrradianceModel
	ssc_data_t data;
	ssc_module_t mod;
	std::unique_ptr<PVIOManager> PVIO;
	std::unique_ptr<IrradianceModel> irradianceModel;

	/// Set up the irradiance case
	void SetUp()
	{
		data = ssc_data_create();
		pvsamv_nofinancial_default(data);

		mod = ssc_module_create(const_cast<char*>("pvsamv2"));
		if (NULL == mod)
		{
			printf("error: could not create 'pvsamv2' module.");
			ssc_data_free(data);
			return;
		}
		if (ssc_module_exec(mod, data) == 0)
		{
			printf("error during simulation.");
			ssc_module_free(mod);
			ssc_data_free(data);
		}
		compute_module *cm = static_cast<compute_module*>(mod);

		// Declare these as smart pointers.  Simply means won't worry about memory deallocation later
		// This line is currently failing, down in core.cpp, due to the cm object not being able find the variables
		std::unique_ptr<PVIOManager> tmp(new PVIOManager(*cm));
		PVIO = std::move(tmp);

		std::unique_ptr<IrradianceModel> tmp2(new IrradianceModel(PVIO.get()));
		irradianceModel = std::move(tmp2);
	}
	/// Tear down case. 
	void TearDown() {
		if (data) {
			ssc_data_clear(data);
		}
	}
};

#endif