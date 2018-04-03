#include <gtest/gtest.h>
#include "lib_pv_irradiance_test.h"

TEST_F(LibPVIrradianceTest, DefaultNoFinancialModel) {

	EXPECT_EQ(irradianceModel->RunSingleStep(0), EXIT_SUCCESS);
}

