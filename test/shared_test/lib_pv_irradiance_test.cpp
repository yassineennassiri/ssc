#include <gtest/gtest.h>
#include "lib_pv_irradiance_test.h"

TEST_F(LibPVIrradianceTest, DefaultNoFinancialModel) {

	// Run with defaults for a step
	EXPECT_EQ(irradianceModel->RunSingleStep(0), EXIT_SUCCESS);

	// Update sky model and run for a step
	PVIO->getIrradianceIO()->skyModel = 0;
	EXPECT_EQ(irradianceModel->RunSingleStep(0), EXIT_SUCCESS);


}

