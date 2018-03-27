#ifndef _LIB_PV_IRRADIANCE_TEST_H_
#define _LIB_PV_IRRADIANCE_TEST_H_


#include <lib_pv_irradiance.h>
#include "../input_cases/pvsamv1_cases.h"


class LibPVIrradiance : public ::testing::Test {

public:

	// Member variables, including PVIOManager and IrradianceModel
	ssc_data_t data;
	ssc_module_t mod;
	PVIOManager * PVIO;
	IrradianceModel * irradianceModel;

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
		compute_module *cm = static_cast<compute_module*>(mod);

		// Declare these as smart pointers.  Simply means won't worry about memory deallocation later
		PVIO = new PVIOManager(*cm);
		irradianceModel = new IrradianceModel(PVIO);
	}
	/// Tear down case. 
	void TearDown() {
		if (data) {
			ssc_data_clear(data);
		}
	}
};

#endif