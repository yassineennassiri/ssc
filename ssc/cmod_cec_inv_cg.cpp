#include "core.h"
#include <vector>
#include <sstream>

#include "lsqfit.h"


  
static var_info vtab_cec_inv_cg[] = {

/*   VARTYPE           DATATYPE         NAME                         LABEL                                           UNITS     META                      GROUP          REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
	{ SSC_INPUT, SSC_NUMBER, "cec_inv_cg_Paco", "Rated max output", "W", "", "", "*", "", "" },

	{ SSC_INPUT, SSC_NUMBER, "cec_inv_cg_power_units", "Sample data units for power output", "0=W,1=kW", "", "", "?=0", "INTEGER,MIN=0,MAX=1", "" },
	// each sample has 18x3 entries:
	// 6 output power percentages 10%, 20%, 30%, 50%, 75%, 100% of rated power
	// 3 voltages Vmin, Vnom, Vmax that the 6 output powers measured at
	// for a total of 18 (=6x3) rows
	// 3 measured values for each row - Output Power, Input Voltage and Efficiency
	{ SSC_INPUT, SSC_MATRIX, "cec_inv_cg_test_samples", "Sample data", "", "", "", "*", "", "" },
	
	/* from pvsamv1
		{ SSC_INPUT,        SSC_NUMBER,      "mppt_low_inverter",                           "Minimum inverter MPPT voltage window",                    "Vdc",     "",                     "pvsamv1",       "",                    "?=0",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "mppt_hi_inverter",                            "Maximum inverter MPPT voltage window",                    "Vdc",     "",                     "pvsamv1",       "",                    "?=0",                              "" },

	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_c0",                                  "Curvature between ac-power and dc-power at ref",          "1/W",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_c1",                                  "Coefficient of Pdco variation with dc input voltage",     "1/V",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_c2",                                  "Coefficient of Pso variation with dc input voltage",      "1/V",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_c3",                                  "Coefficient of Co variation with dc input voltage",       "1/V",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_paco",                                "AC maximum power rating",                                 "Wac",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_pdco",                                "DC input power at which ac-power rating is achieved",     "Wdc",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_pnt",                                 "AC power consumed by inverter at night",                  "Wac",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_pso",                                 "DC power required to enable the inversion process",       "Wdc",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_vdco",                                "DC input voltage for the rated ac-power rating",          "Vdc",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "inv_snl_vdcmax",                              "Maximum dc input operating voltage",                      "Vdc",     "",                     "pvsamv1",       "inverter_model=0",                    "",                              "" },
	*/

	// intermediate outputs for testing and validation
	{ SSC_OUTPUT, SSC_MATRIX, "cec_inv_cg_J87_L91", "Excel cells J87:L91", "", "", "", "*", "", "" },

	// test data reorganized
	{ SSC_OUTPUT, SSC_MATRIX, "cec_inv_cg_Vmin", "Vmin for least squares fit", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_MATRIX, "cec_inv_cg_Vnom", "Vnom for least squares fit", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_MATRIX, "cec_inv_cg_Vmax", "Vmax for least squares fit", "", "", "", "*", "", "" },

	// quadratic fits
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_Vmin_abc", "Vmin a,b,c for least squares fit", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_Vnom_abc", "Vnom a,b,c for least squares fit", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_Vmax_abc", "Vmax a,b,c for least squares fit", "", "", "", "*", "", "" },

	//intermediates based on quadratic least squares
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_Vdc", "Vdc at Vmin, Vnom, Vmax", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_Vdc_Vnom", "Vdc - Vnom at Vmin, Vnom, Vmax", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_Pdco", "Pdco at Vmin, Vnom, Vmax", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_Psco", "Psco at Vmin, Vnom, Vmax", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_C0", "C0 at Vmin, Vnom, Vmax", "", "", "", "*", "", "" },

	// intermediates based on linear least squares
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_C1", "C1 at m and b", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_C2", "C1 at m and b", "", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "cec_inv_cg_C3", "C1 at m and b", "", "", "", "*", "", "" },

	// outputs Pdco, Vdco, Pso, c0, c1, c2, c3
	{ SSC_OUTPUT, SSC_NUMBER, "Pdco", "CEC generated Pdco", "Wac", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "Vdco", "CEC generated Vdco", "Vdc", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "Pso", "CEC generated Pso", "Wdc", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "c0", "CEC generated c0", "1/W", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "c1", "CEC generated c1", "1/V", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "c2", "CEC generated c2", "1/V", "", "", "*", "", "" },
	{ SSC_OUTPUT, SSC_NUMBER, "c3", "CEC generated c3", "1/V", "", "", "*", "", "" },

	var_info_invalid };


double Quadratic_fit_eqn(double _x, double *par, void *)
{
	return par[0] * _x * _x + par[1] * _x + par[2];
}

double Linear_fit_eqn(double _x, double *par, void *)
{
	return par[0] * _x + par[1];
}

class cm_cec_inv_cg : public compute_module
{
private:

public:
	cm_cec_inv_cg()
	{
		add_var_info( vtab_cec_inv_cg);
	}

	void exec( ) throw( general_error )
	{
		size_t count, i, j, nrows, ncols; 
		// rated output ac
		double Paco = as_double("cec_inv_cg_Paco");

		// 6 columns period, tier, max usage, max usage units, buy, sell
		ssc_number_t *cec_inv_cg_test_samples_in = as_matrix("cec_inv_cg_test_samples", &nrows, &ncols);
		if (nrows != 18)
		{
			std::ostringstream ss;
			ss << "The samples table must have 18 rows. Number of rows in samples table provided is " << nrows << " rows.";
			throw exec_error("cec_inv_cg", ss.str());
		}
		if ((ncols % 3) != 0)
		{
			std::ostringstream ss;
			ss << "The samples table must have number of columns divisible by 3. Number of columns in samples table provided is " << ncols << " columns.";
			throw exec_error("cec_inv_cg", ss.str());
		}
		size_t num_samples = ncols / 3;
		size_t columns_per_sample = 3;
		util::matrix_t<float> cec_inv_cg_test_samples(nrows, ncols);
		cec_inv_cg_test_samples.assign(cec_inv_cg_test_samples_in, nrows, ncols);

		util::matrix_t<ssc_number_t> &cec_inv_cg_J87_L91 = allocate_matrix("cec_inv_cg_J87_L91", 5, 3);
		ssc_number_t vdc = 0, Pout = 0, eff = 0, Pin=0, Pin2=0;
		ssc_number_t Vmin_Pout_avg = 0, Vmin_Pin_avg = 0, Vmin_Pin2_avg = 0;
		ssc_number_t Vnom_Pout_avg = 0, Vnom_Pin_avg = 0, Vnom_Pin2_avg = 0;
		ssc_number_t Vmax_Pout_avg = 0, Vmax_Pin_avg = 0, Vmax_Pin2_avg = 0;
		for (size_t k = 0; k < cec_inv_cg_J87_L91.ncols(); k++)
			cec_inv_cg_J87_L91.at(0, k) = 0.0;

		// set Pout, Pin=Pout/eff and Pin^2 for least squares fit
		// 6 is for the required power output percentages at each voltage for each sample
		util::matrix_t<ssc_number_t> &cec_inv_cg_Vmin = allocate_matrix("cec_inv_cg_Vmin", 6 * num_samples, 3);
		util::matrix_t<ssc_number_t> &cec_inv_cg_Vnom = allocate_matrix("cec_inv_cg_Vnom", 6 * num_samples, 3);
		util::matrix_t<ssc_number_t> &cec_inv_cg_Vmax = allocate_matrix("cec_inv_cg_Vmax", 6 * num_samples, 3);


		ssc_number_t *cec_inv_cg_Vdc = allocate("cec_inv_cg_Vdc", 3);
		ssc_number_t *cec_inv_cg_Vdc_Vnom = allocate("cec_inv_cg_Vdc_Vnom", 3);
		ssc_number_t *cec_inv_cg_Pdco = allocate("cec_inv_cg_Pdco", 3);
		ssc_number_t *cec_inv_cg_Psco = allocate("cec_inv_cg_Psco", 3);
		ssc_number_t *cec_inv_cg_C0 = allocate("cec_inv_cg_C0", 3);
		ssc_number_t *cec_inv_cg_C1 = allocate("cec_inv_cg_C1", 2);
		ssc_number_t *cec_inv_cg_C2 = allocate("cec_inv_cg_C2", 2);
		ssc_number_t *cec_inv_cg_C3 = allocate("cec_inv_cg_C3", 2);

		for (i = 0; i < 3; i++)
			cec_inv_cg_Vdc[i] = 0;

		for (j = 0; j < num_samples; j++)
		{
			for (i = 0; i < cec_inv_cg_test_samples.nrows(); i++)
			{
				vdc = cec_inv_cg_test_samples.at(i, j*columns_per_sample + 1);
				Pout = cec_inv_cg_test_samples.at(i, j*columns_per_sample);
				eff = cec_inv_cg_test_samples.at(i, j*columns_per_sample+2);
				Pin = Pout;
				if (eff != 0.0f) Pin = 100.0*Pout / eff;
				Pin2 = Pin*Pin;
				if (i < 6) // Vmin 0 offset
				{
					cec_inv_cg_J87_L91.at(0, 0) += vdc;
					cec_inv_cg_Vdc[0] += vdc;
					cec_inv_cg_Vmin.at(j * 6 + i, 0) = Pout;
					cec_inv_cg_Vmin.at(j * 6 + i, 1) = Pin;
					cec_inv_cg_Vmin.at(j * 6 + i, 2) = Pin2;
					Vmin_Pout_avg += Pout;
					Vmin_Pin_avg += Pin;
					Vmin_Pin2_avg += Pin2;
				}
				else if (i < 12) // Vnom 6 offset 
				{
					cec_inv_cg_J87_L91.at(0, 1) += vdc;
					cec_inv_cg_Vdc[1] += vdc;
					cec_inv_cg_Vnom.at(j * 6 + i - 6, 0) = Pout;
					cec_inv_cg_Vnom.at(j * 6 + i-6, 1) = Pin;
					cec_inv_cg_Vnom.at(j * 6 + i-6, 2) = Pin2;
					Vnom_Pout_avg += Pout;
					Vnom_Pin_avg += Pin;
					Vnom_Pin2_avg += Pin2;
				}
				else // Vmax 12 offset
				{
					cec_inv_cg_J87_L91.at(0, 2) += vdc;
					cec_inv_cg_Vdc[2] += vdc;
					cec_inv_cg_Vmax.at(j * 6 + i - 12, 0) = Pout;
					cec_inv_cg_Vmax.at(j * 6 + i - 12, 1) = Pin;
					cec_inv_cg_Vmax.at(j * 6 + i - 12, 2) = Pin2;
					Vmax_Pout_avg += Pout;
					Vmax_Pin_avg += Pin;
					Vmax_Pin2_avg += Pin2;
				}
			}
		}
		
		


		// do simple least square fit per CoefGenerator Worksheet
		// y=Pout, x1=Pin=Pout/eff, x2=Pin^2
		// 6 output power ranges for each Vmin, Vnom and Vmax for each sample
		Vmin_Pout_avg /= (6 * num_samples);
		Vmin_Pin_avg /= (6 * num_samples);
		Vmin_Pin2_avg /= (6 * num_samples);
		Vnom_Pout_avg /= (6 * num_samples);
		Vnom_Pin_avg /= (6 * num_samples);
		Vnom_Pin2_avg /= (6 * num_samples);
		Vmax_Pout_avg /= (6 * num_samples);
		Vmax_Pin_avg /= (6 * num_samples);
		Vmax_Pin2_avg /= (6 * num_samples);


		ssc_number_t *cec_inv_cg_Vmin_abc = allocate("cec_inv_cg_Vmin_abc", 3);
		ssc_number_t *cec_inv_cg_Vnom_abc = allocate("cec_inv_cg_Vnom_abc", 3);
		ssc_number_t *cec_inv_cg_Vmax_abc = allocate("cec_inv_cg_Vmax_abc", 3);

		/* linear least squares
		ssc_number_t mPin_num = 0, mPin_den = 0, mPin=0;
		ssc_number_t mPin2_num = 0, mPin2_den = 0, mPin2=0;
		ssc_number_t mPout_b = 0;
		for (i = 0; i < cec_inv_cg_Vmin.nrows(); i++)
		{
			// mPin
			mPin_num += (cec_inv_cg_Vmin.at(i, 1) - Vmin_Pin_avg)*(cec_inv_cg_Vmin.at(i, 0) - Vmin_Pout_avg);
			mPin_den += (cec_inv_cg_Vmin.at(i, 1) - Vmin_Pin_avg)*(cec_inv_cg_Vmin.at(i, 1) - Vmin_Pin_avg);
			// mPin2
			mPin2_num += (cec_inv_cg_Vmin.at(i, 2) - Vmin_Pin2_avg)*(cec_inv_cg_Vmin.at(i, 0) - Vmin_Pout_avg);
			mPin2_den += (cec_inv_cg_Vmin.at(i, 2) - Vmin_Pin2_avg)*(cec_inv_cg_Vmin.at(i, 2) - Vmin_Pin2_avg);
		}
		mPin = mPin_num;
		if (mPin_den != 0) mPin /= mPin_den;
		mPin2 = mPin2_num;
		if (mPin2_den != 0) mPin2 /= mPin2_den;
		mPout_b = Vmin_Pout_avg - mPin*Vmin_Pin_avg - mPin2*Vmin_Pin2_avg;

		cec_inv_cg_Vmin_abc[0] = mPin;
		cec_inv_cg_Vmin_abc[1] = mPin2;
		cec_inv_cg_Vmin_abc[2] = mPout_b;
		*/



		std::vector<double> Pout_vec(cec_inv_cg_Vmin.nrows());
		std::vector<double> Pin_vec(cec_inv_cg_Vmin.nrows());
		int info;
		double C[3];// initial guesses for lsqfit

		// Vmin non-linear
		for (i = 0; i < cec_inv_cg_Vmin.nrows(); i++)
		{
			Pin_vec[i] = cec_inv_cg_Vmin.at(i, 1);
			Pout_vec[i] = cec_inv_cg_Vmin.at(i, 0);
		}
		C[0] = -1e-6;
		C[1] = 1;
		C[2] = 1e3;
		if (!(info=lsqfit(Quadratic_fit_eqn, 0, C, 3, &Pin_vec[0], &Pout_vec[0], cec_inv_cg_Vmin.nrows())))
		{
			throw exec_error("cec_inv_cg", util::format("error in nonlinear least squares fit, error %d", info));
			return;
		}
		cec_inv_cg_Vmin_abc[0] = C[0];
		cec_inv_cg_Vmin_abc[1] = C[1];
		cec_inv_cg_Vmin_abc[2] = C[2];

		// Vnom non-linear
		for (i = 0; i < cec_inv_cg_Vnom.nrows(); i++)
		{
			Pin_vec[i] = cec_inv_cg_Vnom.at(i, 1);
			Pout_vec[i] = cec_inv_cg_Vnom.at(i, 0);
		}
		C[0] = -1e-6;
		C[1] = 1;
		C[2] = 1e3;
		if (!(info = lsqfit(Quadratic_fit_eqn, 0, C, 3, &Pin_vec[0], &Pout_vec[0], cec_inv_cg_Vnom.nrows())))
		{
			throw exec_error("cec_inv_cg", util::format("error in nonlinear least squares fit, error %d", info));
			return;
		}
		cec_inv_cg_Vnom_abc[0] = C[0];
		cec_inv_cg_Vnom_abc[1] = C[1];
		cec_inv_cg_Vnom_abc[2] = C[2];

		// Vmax non-linear
		for (i = 0; i < cec_inv_cg_Vmax.nrows(); i++)
		{
			Pin_vec[i] = cec_inv_cg_Vmax.at(i, 1);
			Pout_vec[i] = cec_inv_cg_Vmax.at(i, 0);
		}
		C[0] = -1e-6;
		C[1] = 1;
		C[2] = 1e3;
		if (!(info = lsqfit(Quadratic_fit_eqn, 0, C, 3, &Pin_vec[0], &Pout_vec[0], cec_inv_cg_Vmax.nrows())))
		{
			throw exec_error("cec_inv_cg", util::format("error in nonlinear least squares fit, error %d", info));
			return;
		}
		cec_inv_cg_Vmax_abc[0] = C[0];
		cec_inv_cg_Vmax_abc[1] = C[1];
		cec_inv_cg_Vmax_abc[2] = C[2];

		// Fill in intermediate values
		//Vdc (Vmin, Vnom, Vmax)
		for (i = 0; i < 3;i++)
			cec_inv_cg_Vdc[i] /= (6 * num_samples);
		for (size_t k = 0; k < cec_inv_cg_J87_L91.ncols(); k++)
			cec_inv_cg_J87_L91.at(0, k) /= (6 * num_samples);

		// Vdc-Vnom
		for (i = 0; i < 3; i++)
			cec_inv_cg_Vdc_Vnom[i] = cec_inv_cg_Vdc[i] - cec_inv_cg_Vdc[1];
		for (size_t k = 0; k < cec_inv_cg_J87_L91.ncols(); k++)
			cec_inv_cg_J87_L91.at(1, k) = cec_inv_cg_J87_L91.at(0, k) - cec_inv_cg_J87_L91.at(0, 1);

		ssc_number_t a, b, c;
		// Pdco and Psco and C0
		a = cec_inv_cg_Vmin_abc[0];
		b = cec_inv_cg_Vmin_abc[1];
		c = cec_inv_cg_Vmin_abc[2];
		cec_inv_cg_Pdco[0] = (-b + sqrt(b*b - 4 * a*(c - Paco)));
		cec_inv_cg_Psco[0] = (-b + sqrt(b*b - 4 * a*c));
		cec_inv_cg_C0[0] = a;
		if (a != 0)
		{
			cec_inv_cg_Pdco[0] /= (2.0*a);
			cec_inv_cg_Psco[0] /= (2.0*a);
		}

		a = cec_inv_cg_Vnom_abc[0];
		b = cec_inv_cg_Vnom_abc[1];
		c = cec_inv_cg_Vnom_abc[2];
		cec_inv_cg_Pdco[1] = (-b + sqrt(b*b - 4 * a*(c - Paco)));
		cec_inv_cg_Psco[1] = (-b + sqrt(b*b - 4 * a*c));
		cec_inv_cg_C0[1] = a;
		if (a != 0)
		{
			cec_inv_cg_Pdco[1] /= (2.0*a);
			cec_inv_cg_Psco[1] /= (2.0*a);
		}

		// TODO - limit Psco max to not be less than zero per note in Workbook
		a = cec_inv_cg_Vmax_abc[0];
		b = cec_inv_cg_Vmax_abc[1];
		c = cec_inv_cg_Vmax_abc[2];
		cec_inv_cg_Pdco[2] = (-b + sqrt(b*b - 4 * a*(c - Paco)));
		cec_inv_cg_Psco[2] = (-b + sqrt(b*b - 4 * a*c));
		cec_inv_cg_C0[2] = a;
		if (a != 0)
		{
			cec_inv_cg_Pdco[2] /= (2.0*a);
			cec_inv_cg_Psco[2] /= (2.0*a);
		}

		// C1, C2, C3 linear least squares
		// C1 Y=Pdco, X=Vdc-Vnom
		std::vector<double> X(3);
		std::vector<double> Y(3);
		double slope, intercept;

		// C1 using linear least squares fit
		for (i = 0; i < 3; i++)
		{
			X[i] = cec_inv_cg_Vdc_Vnom[i];
			Y[i] = cec_inv_cg_Pdco[i];
		}
		if ((info = linlsqfit(&slope, &intercept, &X[0], &Y[0], 3)))
		{
			throw exec_error("cec_inv_cg", util::format("error in linear least squares fit, error %d", info));
			return;
		}
		cec_inv_cg_C1[0] = slope;
		cec_inv_cg_C1[1] = intercept;


		// C2 using linear least squares fit
		for (i = 0; i < 3; i++)
		{
			X[i] = cec_inv_cg_Vdc_Vnom[i];
			Y[i] = cec_inv_cg_Psco[i];
		}
		if ((info = linlsqfit(&slope, &intercept, &X[0], &Y[0], 3)))
		{
			throw exec_error("cec_inv_cg", util::format("error in linear least squares fit, error %d", info));
			return;
		}
		cec_inv_cg_C2[0] = slope;
		cec_inv_cg_C2[1] = intercept;

		// C2 using linear least squares fit
		for (i = 0; i < 3; i++)
		{
			X[i] = cec_inv_cg_Vdc_Vnom[i];
			Y[i] = cec_inv_cg_C0[i];
		}
		if ((info = linlsqfit(&slope, &intercept, &X[0], &Y[0], 3)))
		{
			throw exec_error("cec_inv_cg", util::format("error in linear least squares fit, error %d", info));
			return;
		}
		cec_inv_cg_C3[0] = slope;
		cec_inv_cg_C3[1] = intercept;

		// vdco is the average of Vnom of all samples column 2 and rows 7 through 12

		assign("Pdco", (var_data)cec_inv_cg_C1[1]);
		assign("Vdco", (var_data)cec_inv_cg_Vdc[1]);
		assign("Pso", (var_data)cec_inv_cg_C2[1]);
		assign("c0", (var_data)cec_inv_cg_C3[1]);
		assign("c1", (var_data)(cec_inv_cg_C1[0] / cec_inv_cg_C1[1]));
		assign("c2", (var_data)(cec_inv_cg_C2[0] / cec_inv_cg_C2[1]));
		assign("c3", (var_data)(cec_inv_cg_C3[0] / cec_inv_cg_C3[1]));
	}





};

DEFINE_MODULE_ENTRY( cec_inv_cg, "CEC Inverter Coefficient Generator", 1 );


