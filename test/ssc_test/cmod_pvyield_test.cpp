#include <gtest/gtest.h>

#include "cmod_pvyield_test.h"
#include "../input_cases/pvyield_cases.h"
#include "../input_cases/weather_inputs.h"


//using std::cout;
//using std::endl;


/// Test PVSAMv1 with inputs from PVYield
TEST_F(PVYieldTimo, DefaultTimoModel)
{

	ssc_data_t data = ssc_data_create();
	int pvsam_errors = pvyield_test(data);
	EXPECT_FALSE(pvsam_errors);
	printf("ssc version %d build information %s", ssc_version(), ssc_build_info());

	if (!pvsam_errors)
	{
		ssc_number_t annual_energy;
		ssc_data_get_number(data, "annual_energy", &annual_energy);
		EXPECT_NEAR(annual_energy, 7407907, 7407907e-4) << "Annual energy.";

		ssc_number_t capacity_factor;
		ssc_data_get_number(data, "capacity_factor", &capacity_factor);
		EXPECT_NEAR(capacity_factor, 20.219496, m_error_tolerance_lo) << "Capacity factor";

		ssc_number_t kwh_per_kw;
		ssc_data_get_number(data, "kwh_per_kw", &kwh_per_kw);
		EXPECT_NEAR(kwh_per_kw, 1771.227905, m_error_tolerance_hi) << "Energy yield";

		ssc_number_t performance_ratio;
		ssc_data_get_number(data, "performance_ratio", &performance_ratio);
		EXPECT_NEAR(performance_ratio, -14.485646, m_error_tolerance_lo) << "Energy yield";

	/*
	cout << "-----------------------------------------------------" << endl;
	cout << "annual_energy (7407907):                 " << annual_energy << endl;
	cout << "capacity_factor (20.219496):          " << capacity_factor << endl;
	cout << "kwh_per_kw (1771.227905):           " << kwh_per_kw << endl;
	cout << "performance_ratio (-14.485646):                 " << performance_ratio << endl;
	cout << "-----------------------------------------------------" << endl;
	*/
	}

}
