#include <gtest/gtest.h>
#include "lib_pv_irradiance_test.h"

TEST_F(LibPVIrradiance, DefaultNoFinancialModel) {

	irradianceModel->RunSingleStep(0);
}

