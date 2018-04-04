#ifndef _LIB_PV_IRRADIANCE_TEST_H_
#define _LIB_PV_IRRADIANCE_TEST_H_


#include <lib_pv_irradiance.h>
#include <cmod_pvsamv2.h>
#include "lib_pv_io.h"
#include "core.h"

#include "../input_cases/pvsamv1_cases.h"
#include "lib_pv_io_test.h"

class fakeComputeModule : public compute_module {
	ssc_data_t data;
public:
	fakeComputeModule() {
		data = ssc_data_create();
		/// temporarily use report generator data
		pvsamv_nofinancial_default(data);
		var_table *vt = static_cast<var_table*>(data);
		m_vartab = vt;
	}
	void exec() {}
};

class LibPVIrradianceTest : public ::testing::Test{
public:
	fakeComputeModule* cm;
	IrradianceModel* irradianceModel;
	PVIOManager* PVIO;
	void SetUp() {
		cm = new fakeComputeModule();
		PVIO = new PVIOManager(cm);
		irradianceModel = new IrradianceModel(PVIO);
	}
	void TearDown() {
		delete PVIO;
		delete irradianceModel;
	}
};

#endif