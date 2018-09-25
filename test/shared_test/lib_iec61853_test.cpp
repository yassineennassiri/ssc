#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "lib_iec61853_test.h"
#include <lib_irradproc.h>


TEST_F(IEC61853Test, InitialGuessesFrom61215) {
	std::vector<std::vector<double>> params6parvs5Par = { {1.10871,7.10E-13,10.3568,264.923,2.14945, 1.09209,1.47E-08,7.72842,307.62,3.33896},
		{1.18283, 4.35E-13, 10.0181, 237.744, 2.13429, 1.16457, 9.60E-09, 7.5535, 272.796, 3.2926},
		{4.81852, 1.83E-11, 1.13416, 22.1686, 0.883573 , 4.71262, 9.00E-07, 0.841604, 27.8294, 1.50673},
		{4.63934, 2.52E-11, 1.10212, 23.4403, 0.902403 , 4.55162, 4.16E-07, 0.830463, 28.283, 1.44845},
		{1.23955, 6.58E-15, 15.8526, 445.972, 2.68498, 1.22599, 3.66E-11, 13.8399, 505.96, 3.64539},
		{1.23692, 6.28E-15, 16.1935, 429.746, 2.63714, 1.22275, 3.83E-11, 14.1845, 489.786, 3.59246},
		{2.53002, 8.37E-12, 1.95798, 196.028, 1.6052, 2.51614, 3.11E-08, 1.32404, 274.564, 2.33069},
		{6.41606, 3.20E-10, 1.7364, 23.1282, 1.73843, 6.13996, 0.00014549, 1.0163, 33.7909, 3.89953},
		{6.6911, 4.98E-10, 1.63197, 25.8665, 1.79502, 6.44946, 9.01E-05, 1.03099, 39.2835, 3.76087},
		{2.57861, 6.46E-12, 2.2845, 118.89, 1.58976, 2.55333, 9.57E-08, 1.48713, 150.951, 2.48656},
		{5.59693, 5.66E-12, 0.586141, 253.221, 1.84825, 5.59379, 4.24E-11, 0.525622, 281.032, 1.99359},
		{5.5408, 8.32E-12, 0.591446, 371.726, 1.84598, 5.53771, 1.08E-10, 0.51189, 468.647, 2.03747},
		{2.74615, 6.43E-11, 0.427333, 227.482, 0.902943 , 2.74256, 2.69E-09, 0.296496, 313.895, 1.06531},
		{2.75661, 5.53E-11, 0.479298, 199.524, 0.897479 , 2.75197, 1.59E-09, 0.365917, 251.726, 1.03902},
		{2.74581, 5.66E-11, 0.478258, 225.685, 0.896246 , 2.74137, 2.73E-09, 0.344366, 316.114, 1.06367},
		{2.74589, 6.08E-11, 0.507344, 235.96, 0.898401 , 2.7407, 2.15E-09, 0.386535, 329.595, 1.05097},
		{5.07577, 1.38E-10, 0.369558, 158.954, 0.891735 , 5.05389, 1.58E-08, 0.281059, 347.676, 1.10732},
		{5.11485, 1.08E-10, 0.400293, 121.093, 0.883824 , 5.10101, 6.45E-09, 0.325849, 186.931, 1.0598},
		{5.13103, 8.15E-12, 0.528116, 46.9893, 0.813009 , 5.10682, 4.32E-09, 0.416537, 57.5364, 1.05757},
		{5.13834, 1.13E-10, 0.37723, 86.4036, 0.900452 , 5.12876, 1.73E-09, 0.327821, 101.194, 1.01296} };

	iec61853_module_t solver;
	for (size_t i = 0; i < mmSTCVector.size(); i++) {
		moduleMeasurements mm = mmSTCVector[i];
		/*initial guesses */
		double Il = params6parvs5Par[i][0];
		double Io = params6parvs5Par[i][1];
		double Rs = params6parvs5Par[i][2];
		double Rsh = params6parvs5Par[i][3];
		double a = params6parvs5Par[i][4];

		bool solved = solver.solve(mm.Voc, mm.Isc, mm.Vmp, mm.Imp, a, &Il, &Io, &Rs, &Rsh);
		double RMSError = pow( 
			pow(Il - params6parvs5Par[i][0],2) +
			pow(Io - params6parvs5Par[i][1], 2) +
			pow(Rs - params6parvs5Par[i][2], 2) +
			pow(Rsh - params6parvs5Par[i][3], 2)
			,.5);

		printf("module: %s\tsolved? %d \t error: %f \n", mm.name.c_str(), solved, RMSError);
	}
}



TEST_F(IEC61853Test, ParameterEstimateWithIECSolverTest) {
	for (size_t m = 0; m < 20; m++) {
		std::vector<double> testData = groupByModule(m);
		util::matrix_t<double> testMatrix(18, 6, &testData);
		util::matrix_t<double> par;

		if (!solver.modIEC.calculate(testMatrix, mmVector[m*18].nSer, mmVector[m*18].type, par, true)) continue;

		for (size_t n = 0; n < 18; n++) {
			mmVector[m * 18 + n].Il = par.at(n, 0);
			mmVector[m * 18 + n].Io = par.at(n, 1);
			mmVector[m * 18 + n].Rs = par.at(n, 2);
			mmVector[m * 18 + n].Rsh = par.at(n, 3);
			mmVector[m * 18 + n].A = par.at(n, 4);
		}

	}

	for (size_t i = 0; i < mmVector.size(); i++) {
		moduleMeasurements* mm = &mmVector[i];
		int tech_id = mm->type;
		double Vmp = mm->Vmp;
		double Imp = mm->Imp;
		double Voc = mm->Voc;
		double Isc = mm->Isc;
		double alpha = 0;
		double beta = 0;
		double gamma = 0;
		int nser = mmVector[i].nSer;

		module6par m(tech_id, Vmp, Imp, Voc, Isc, beta, alpha, gamma, nser, 298.15);
		m.setUncertaintyPmp(5);
		//std::cout << "module " << i << ": ";
	
		m.Io = mm->Io;
		m.Il = mm->Il;
		m.Rs = mm->Rs;
		m.Rsh = mm->Rsh;
		m.a = mm->A;
		m.Adj = 0;

		mm->solved = m.sanity();

		double vStep;
		if (i % 18 == 0) {
			vStep = mmVector[i].Voc / 10.;
			for (size_t j = 0; j < 17; j++) vStep = std::max(vStep, mmVector[i + j].Voc / 10);
		}
		else {
			size_t ref = i - (size_t)(i % 18);
			vStep = mmVector[ref].currentAt10Voltages[1][0];
		}

		for (size_t s = 0; s < 10; s++) {
			double v = vStep * (double)s;
			std::vector<double> iv;
			iv.push_back(v);
			iv.push_back(module6par::current(v, m.Il, m.Io, m.Rs, m.a, m.Rsh, mm->Imp));
			mm->currentAt10Voltages.push_back(iv);
		}

	}

	mmVectorToCSV();
}

TEST_F(IEC61853Test, CheckTempCoeffs) {
	std::string outputFileName = "C:/Users/dguittet/Documents/IEC 61853 Modeling/Data For Validating Models/tempCoeffsIEC61853.csv";
	std::ofstream outputFile;
	outputFile.open(outputFileName);
	EXPECT_TRUE(outputFile.is_open());

	for (size_t m = 0; m < 20; m++) {
		std::vector<double> testData = groupByModule(m);
		util::matrix_t<double> testMatrix(18, 6, &testData);
		util::matrix_t<double> par;

		double Pmp0;
		for (size_t i = 0; i < testMatrix.nrows(); i++) {
			if (testMatrix(i, iec61853_module_t::COL_TC) == 25. && testMatrix(i, iec61853_module_t::COL_IRR) == 1000.) {
				Pmp0 = testMatrix(i, iec61853_module_t::COL_PMP);
			}
		}

		double betaVoc, alphaIsc, gammaPmp;
		// estimate beta VOC at STC irradiance (1000 W/m2)
		EXPECT_FALSE(!solver.modIEC.tcoeff(testMatrix, 4, 1000.0, &betaVoc, false));
			

		// estimate the alpha ISC at STC conditions
		EXPECT_FALSE(!solver.modIEC.tcoeff(testMatrix, 5, 1000.0, &alphaIsc, false));

		// estimate the gamma PMP at STC conditions
		EXPECT_FALSE(!solver.modIEC.tcoeff(testMatrix, 2, 1000.0, &gammaPmp, false));

		gammaPmp *= 100.0 / (Pmp0);

		outputFile << alphaIsc << ", " << betaVoc << ", " << gammaPmp << "\n";

	}
	outputFile.close();
}

TEST_F(IEC61215Test, unsolvedModules) {
	// Uncertainty for crystalline silicon modules: 
	// Pm = ± 2.8%,  Isc = ± 2.3%, Imp = ± 2.3%, Voc = ± 0.3%, Vmp = ± 0.7%
	for (size_t i = 0; i < mmVector.size(); i++) {
		if (mmVector[i].solved == 0) continue;

		double Io_best = mmVector[i].Io;
		double Il_best = mmVector[i].Il;
		double Rs_best = mmVector[i].Rs;
		double Rsh_best = mmVector[i].Rsh;
		double Adj_best = mmVector[i].Adj;
		double a_best = mmVector[i].A;
		double solved_best = mmVector[i].solved;
		double pmpCalc = mmVector[i].Vmp* module6par::current(mmVector[i].Vmp, Il_best,
			Io_best, Rs_best, a_best, Rsh_best, mmVector[i].Imp);
		double pmpError = (pmpCalc - mmVector[i].Pm)/mmVector[i].Pm;
		double iocError = module6par::current(mmVector[i].Voc, Il_best,
			Io_best, Rs_best, a_best, Rsh_best, mmVector[i].Imp);
		double bestError = abs(pmpError) + abs(iocError);

		int intervals = 4;
		double uncertainty = 0.030;
		double dPm = mmVector[i].Pm * uncertainty / (double)(intervals/2.);
		double dIsc = mmVector[i].Isc * uncertainty / (double)(intervals / 2.);
		double dImp = mmVector[i].Imp * uncertainty / (double)(intervals / 2.);
		double dVoc = mmVector[i].Voc * uncertainty / (double)(intervals/2);
		double dVmp = mmVector[i].Vmp * uncertainty / (double)(intervals/2);
		double Pm0 = mmVector[i].Pm;
		double Isc0 = mmVector[i].Isc;
		double Voc0 = mmVector[i].Voc;
		double Vmp0 = mmVector[i].Vmp;
		double Imp0 = mmVector[i].Imp;
		mmVector[i].Isc -= mmVector[i].Isc * uncertainty;
		mmVector[i].Imp -= mmVector[i].Imp * uncertainty;
		mmVector[i].Voc -= mmVector[i].Voc * uncertainty;
		mmVector[i].Vmp -= mmVector[i].Vmp * uncertainty;
		for (size_t q = 0; q < intervals + 1; q++) {
			if (q>0) mmVector[i].Isc += dIsc;
			for (size_t r = 0; r < intervals + 1; r++) {
				if (r>0) mmVector[i].Imp += dImp;
				if (mmVector[i].Imp > mmVector[i].Isc) continue;
				for (size_t s = 0; s < intervals+1; s++) {
					if (s>0) mmVector[i].Voc += dVoc;
					for (size_t t = 0; t < intervals + 1; t++) {
						if (t>0 )mmVector[i].Vmp += dVmp;
						mmVector[i].Pm = mmVector[i].Vmp * mmVector[i].Imp;
						if (abs((mmVector[i].Pm - Pm0) / Pm0) > uncertainty) continue;
						if (mmVector[i].Vmp > mmVector[i].Voc) continue;
						int tech_id = 2;
						double Vmp = mmVector[i].Vmp;
						double Imp = mmVector[i].Imp;
						double Voc = mmVector[i].Voc;
						double Isc = mmVector[i].Isc;
						int nser = mmVector[i].nSer;
						double alpha = mmVector[i].alphaIsc;
						double beta = mmVector[i].betaVoc;
						double gamma = mmVector[i].gammaPmp;

						// monoSi, multiSi, CdTe, CIS, CIGS, Amorphous 
						module6par m(tech_id, Vmp, Imp, Voc, Isc, beta, alpha, gamma, nser, 298.15);
						m.setUncertaintyPmp(uncertainty*100);
						int attempts = 0;
						double tol = 1e-8;
						std::cout << mmVector[i].Imp << ", " << mmVector[i].Vmp << ", ";
						std::cout << mmVector[i].Isc << ", " << mmVector[i].Voc << "," << std::endl;
						while (attempts < 4 && mmVector[i].solved != 0) {
							mmVector[i].solved = m.solve_with_sanity_and_heuristics<double>(150, tol);
							tol *= 10;
							attempts++;

							mmVector[i].Il = m.Il;
							mmVector[i].Io = m.Io;
							mmVector[i].Rs = m.Rs;
							mmVector[i].Rsh = m.Rsh;
							mmVector[i].A = m.a;
							mmVector[i].Adj = m.Adj;


							double pmpCal = mmVector[i].Vmp* module6par::current(mmVector[i].Vmp, Il_best,
								Io_best, Rs_best, a_best, Rsh_best, mmVector[i].Imp);
							double pmpErr = (pmpCal - mmVector[i].Pm) / mmVector[i].Pm;
							double iocErr = module6par::current(mmVector[i].Voc, Il_best,
								Io_best, Rs_best, a_best, Rsh_best, mmVector[i].Imp);
							double error = abs(pmpErr) + abs(iocErr);

							if (error < bestError) {
								Io_best = mmVector[i].Io;
								Il_best = mmVector[i].Il;
								Rs_best = mmVector[i].Rs;
								Rsh_best = mmVector[i].Rsh;
								Adj_best = mmVector[i].Adj;
								a_best = mmVector[i].A;
								solved_best = mmVector[i].solved;
								bestError = error;
							}
							if (mmVector[i].solved == 0) break;
							std::cout << ", " << q <<", " << r << ", " << s << ", " << t << "\n";
						}
					
						if (mmVector[i].solved == 0) break;	
					}
					if (mmVector[i].solved == 0) break;
					mmVector[i].Vmp = Vmp0 - Vmp0 * uncertainty;
					
				}
				if (mmVector[i].solved == 0) break;
				mmVector[i].Voc = Voc0 - Voc0 * uncertainty;
			}
			mmVector[i].Imp = Imp0 - Imp0 * uncertainty;
			if (mmVector[i].solved == 0) break;
		}
		if (mmVector[i].solved != 0) {
			mmVector[i].Io = Io_best;
			mmVector[i].Il = Il_best;
			mmVector[i].Rs = Rs_best;
			mmVector[i].Rsh = Rsh_best;
			mmVector[i].Adj = Adj_best;
			mmVector[i].A = a_best;
			mmVector[i].solved = (int)solved_best;
		}
	}
	mmVectorToCSV();
}

TEST_F(IEC61215Test, IVCurves) {
	std::ofstream file;
	file.open("C:/Users/dguittet/Documents/IEC 61853 Modeling/Data For Validating Models/moduleIVCurves.csv");
	EXPECT_TRUE(file.is_open());

	for (size_t i = 0; i < mmVector.size(); i++) {
		if (mmVector[i].solved == 0) continue;

		file << "Name," << mmVector[i].name << ", Imp," << mmVector[i].Imp << ", Vmp,";
		file << mmVector[i].Vmp << ", Pmp" << mmVector[i].gammaPmp << ", Voc,";
		file << mmVector[i].Voc << ", Imp, " << mmVector[i].Imp << "\n";

		double vStep = mmVector[i].Voc / 100.;

		std::vector<double> V;
		std::vector<double> I;
		for (size_t s = 0; s < 100; s++) {
			double v = vStep * (double)s;
			V.push_back(v);
			I.push_back(module6par::current(v, mmVector[i].Il, mmVector[i].Io, mmVector[i].Rs, 
				mmVector[i].A, mmVector[i].Rsh, mmVector[i].Imp));
		}
		for (size_t s = 0; s < 100; s++) {
			file << V[s] << ",";
		}
		file << "\n";
		for (size_t s = 0; s < 100; s++) {
			file << I[s] << ",";
		}
		file << "\n\n";
	}
	file.close();
}

TEST_F(IEC61215Test, solveCoefs) {
	for (size_t i = 0; i < mmVector.size(); i++) {
		int tech_id = 2;
		double Vmp = mmVector[i].Vmp;
		double Imp = mmVector[i].Imp;
		double Voc = mmVector[i].Voc;
		double Isc = mmVector[i].Isc;
		int nser = mmVector[i].nSer;
		double alpha = mmVector[i].alphaIsc;
		double beta = mmVector[i].betaVoc;
		double gamma = mmVector[i].gammaPmp;

		// monoSi, multiSi, CdTe, CIS, CIGS, Amorphous 
		module6par m(tech_id, Vmp, Imp, Voc, Isc, beta, alpha, gamma, nser, 298.15);
		m.setUncertaintyPmp(5);
		int attempts = 0;
		double tol = 1e-8;
		while (attempts < 10 && mmVector[i].solved != 0) {
			mmVector[i].solved = m.solve_with_sanity_and_heuristics<double>(300, tol);
			tol *= 10;
			attempts++;
		}

		mmVector[i].Il = m.Il;
		mmVector[i].Io = m.Io;
		mmVector[i].Rs = m.Rs;
		mmVector[i].Rsh = m.Rsh;
		mmVector[i].A = m.a;
		mmVector[i].Adj = m.Adj;
	}
	mmVectorToCSV();
}

TEST_F(IEC61853Test, testMeasurementsAbsorbedIrradiance) {
	// cocoa, eugene, golden
	double tilt[3] = { 28.5, 44, 40};
	double lon[3] = { -80.46, -123.07, -105.18};
	double lat[3] = { 28.39, 44.05, 39.74};
	double tz[3] = { -5, -8, -7 };
	for (size_t m = 0; m < 20; m++) {
		std::string modulename[20] = { "aSiTandem72-46", "aSiTandem90-31", "aSiTriple28324", "aSiTriple28325", "CdTe75638", "CdTe75669",
			"CIGS1-001", "CIGS8-001", "CIGS39013", "CIGS39017", "HIT05662", "HIT05667", "mSi0166", "mSi0188", "mSi0247",
			"mSi0251", "mSi460A8", "mSi460BB", "xSi11246", "xSi12922" };
		std::string filename = "C:\\Users\\dguittet\\Documents\\IEC 61853 Modeling\\Data For Validating Models\\ResultsComparingModels\\"
			+ modulename[m] + ".csv";
		std::vector<std::string> fileLines = getTestMeasurementData(filename);

		std::string outputFileName = filename;
		std::ofstream outputFile;
		outputFile.open(outputFileName);
		EXPECT_TRUE(outputFile.is_open());
		outputFile << fileLines[0] << "POA_calc";

		for (size_t i = 0; i < testMeasurements.size(); i++) {
			std::string name = testMeasurements[i][0].c_str();
			std::string loc = testMeasurements[i][1].c_str();
			double month = atof(testMeasurements[i][2].c_str());
			double day = atof(testMeasurements[i][3].c_str());
			double hour = atof(testMeasurements[i][4].c_str());
			double minute = atof(testMeasurements[i][5].c_str());
			double POA_meas = atof(testMeasurements[i][6].c_str());
			double T = atof(testMeasurements[i][7].c_str());
			double dn = atof(testMeasurements[i][11].c_str());
			double gh = atof(testMeasurements[i][12].c_str());
			double dh = atof(testMeasurements[i][13].c_str());
			size_t locIndex = 2;
			if (loc == "Cocoa") locIndex = 0;
			else if (loc == "Eugene") locIndex = 1;

			irrad irr;
			irr.set_time(2011, (int)month, (int)day, (int)hour, minute, IRRADPROC_NO_INTERPOLATE_SUNRISE_SUNSET);
			irr.set_location(lat[locIndex], lon[locIndex], tz[locIndex]);

			irr.set_sky_model(1, 0.2);
			irr.set_beam_diffuse(dn, dh);
		
			EXPECT_FALSE(irr.calc());
		}

	}
}

TEST_F(IEC61853Test, testMeasurements) {
	// cocoa, eugene, golden
	double tilt[3] = { 28.5, 44, 40 };
	double lon[3] = { -80.46, -123.07, -105.18 };
	double lat[3] = { 28.39, 44.05, 39.74 };
	double tz[3] = { -5, -8, -7 };
	double elev[3] = { 12, 145, 1798};
	for (size_t m = 0; m < 20; m++) {
		std::string modulename[20] = { "aSiTandem72-46", "aSiTandem90-31", "aSiTriple28324", "aSiTriple28325", "CdTe75638", "CdTe75669",
			"CIGS1-001", "CIGS8-001", "CIGS39013", "CIGS39017", "HIT05662", "HIT05667", "mSi0166", "mSi0188", "mSi0247",
			"mSi0251", "mSi460A8", "mSi460BB", "xSi11246", "xSi12922" };
		std::string filename = "C:\\Users\\dguittet\\Documents\\IEC 61853 Modeling\\Data For Validating Models\\ResultsComparingModels\\"
			+ modulename[m] + ".csv";
		std::vector<std::string> fileLines = getTestMeasurementData(filename);

		std::string outputFileName = filename;
		std::ofstream outputFile;
		outputFile.open(outputFileName);
		EXPECT_TRUE(outputFile.is_open());

		util::matrix_t<double> input(18,6);
		for (int c = 0; c < 18; c++) {
			input.set_value(mmVector[m * 18 + c].irradiance, c, 0);
			input.set_value(mmVector[m * 18 + c].temp, c, 1);
			input.set_value(mmVector[m * 18 + c].Pm, c, 2);
			input.set_value(mmVector[m * 18 + c].Vmp, c, 3);
			input.set_value(mmVector[m * 18 + c].Voc, c, 4);
			input.set_value(mmVector[m * 18 + c].Isc, c, 5);
		}
		int nser = mmVector[m * 18].nSer;
		int type = mmVector[m * 18].type;
		util::matrix_t<double> par;

		iec61853_module_t solver;
		solver.AMA[0] = 0.9417;
		solver.AMA[1] = 0.06516;
		solver.AMA[2] = -0.02022;
		solver.AMA[3] = 0.00219;
		solver.AMA[4] = -9.1e-05;
		solver.GlassAR = 0;
		if (!solver.calculate(input, nser, type, par, false)) continue;

		outputFile << fileLines[0] << ",11parIrr,11parPOA\n";
		for (size_t i = 0; i < testMeasurements.size(); i++) {
			std::string name = testMeasurements[i][0].c_str();
			std::string loc = testMeasurements[i][1].c_str();
			double month = atof(testMeasurements[i][2].c_str());
			double day = atof(testMeasurements[i][3].c_str());
			double hour = atof(testMeasurements[i][4].c_str());
			double minute = atof(testMeasurements[i][5].c_str());
			double POA_meas = atof(testMeasurements[i][6].c_str());
			double T = atof(testMeasurements[i][7].c_str());	// pv back-surface temp
			double dn = atof(testMeasurements[i][11].c_str());
			double gh = atof(testMeasurements[i][12].c_str());
			double dh = atof(testMeasurements[i][13].c_str());
			size_t locIndex = 2;
			if (loc == "Cocoa") locIndex = 0;
			else if (loc == "Eugene") locIndex = 1;

			irrad irr;
			irr.set_time(2011, (int)month, (int)day, (int)hour, minute, IRRADPROC_NO_INTERPOLATE_SUNRISE_SUNSET);
			irr.set_location(lat[locIndex], lon[locIndex], tz[locIndex]);

			irr.set_sky_model(1, 0.2);
			irr.set_beam_diffuse(dn, dh);

			irr.set_surface(0, tilt[locIndex], 180., 0., 0, 0.3);

			EXPECT_FALSE(irr.calc());

			// using HDKR model
			int sunup = 0;
			double solazi = 0, solzen = 0, solalt = 0;
			double ibeam, iskydiff, ignddiff;
			double aoi, stilt, sazi, rot, btd;

			irr.get_sun(&solazi, &solzen, &solalt, 0, 0, 0, &sunup, 0, 0, 0);
			irr.get_angles(&aoi, &stilt, &sazi, &rot, &btd);
			irr.get_poa(&ibeam, &iskydiff, &ignddiff, 0, 0, 0);


			pvinput_t in(ibeam, iskydiff, ignddiff, 0, 0,
				T, 0.0, 0.0, 0.0, 0.0,
				solzen, aoi, elev[locIndex],
				stilt, sazi,
				((double)hour) + minute / 60.0,
				Irradiance_IO::DN_DF, false);

			pvoutput_t output;
			solver(in, T, -1, output);
			outputFile << fileLines[i + 1] << "," << output.Power << ",";

			// using measured POA
			pvinput_t in2(0., 0., 0., 0, POA_meas,
				T, 0.0, 0.0, 0.0, 0.0,
				solzen, aoi, elev[locIndex],
				stilt, sazi,
				((double)hour) + minute / 60.0,
				3, true);
			solver(in2, T, -1, output);
			outputFile <<  output.Power << "\n";
		}
		outputFile.close();
	}
}

TEST_F(IEC61215Test, testMeasurements) {
	// cocoa, eugene, golden
	double tilt[3] = { 28.5, 44, 40 };
	double lon[3] = { -80.46, -123.07, -105.18 };
	double lat[3] = { 28.39, 44.05, 39.74 };
	double tz[3] = { -5, -8, -7 };
	double elev[3] = { 12, 145, 1798 };
	for (size_t m = 0; m < 20; m++) {
		std::string modulename[20] = { "aSiTandem72-46", "aSiTandem90-31", "aSiTriple28324", "aSiTriple28325", "CdTe75638", "CdTe75669",
			"CIGS1-001", "CIGS8-001", "CIGS39013", "CIGS39017", "HIT05662", "HIT05667", "mSi0166", "mSi0188", "mSi0247",
			"mSi0251", "mSi460A8", "mSi460BB", "xSi11246", "xSi12922" };
		std::string filename = "C:\\Users\\dguittet\\Documents\\IEC 61853 Modeling\\Data For Validating Models\\ResultsComparingModels\\" 
			+ modulename[m] +".csv";
		std::vector<std::string> fileLines = getTestMeasurementData(filename);

		std::string outputFileName = filename;
		std::ofstream outputFile;
		outputFile.open(outputFileName);
		EXPECT_TRUE(outputFile.is_open());
		outputFile << fileLines[0] << ",6parIrr, 6parPOA\n";

		moduleMeasurements currentModule;
		std::vector<double> powerPredicted;
		for (size_t i = 0; i < testMeasurements.size(); i++) {
			std::string name = testMeasurements[i][0].c_str();
			std::string loc = testMeasurements[i][1].c_str();
			double month = atof(testMeasurements[i][2].c_str());
			double day = atof(testMeasurements[i][3].c_str());
			double hour = atof(testMeasurements[i][4].c_str());
			double minute = atof(testMeasurements[i][5].c_str());
			double POA_meas = atof(testMeasurements[i][6].c_str());
			double T_cell = atof(testMeasurements[i][7].c_str());	// pv back-surface temp
			double dn = atof(testMeasurements[i][11].c_str());
			double gh = atof(testMeasurements[i][12].c_str());
			double dh = atof(testMeasurements[i][13].c_str());
			size_t locIndex = 2;
			if (loc == "Cocoa") locIndex = 0;
			else if (loc == "Eugene") locIndex = 1;

			irrad irr;
			irr.set_time(2011, (int)month, (int)day, (int)hour, minute, IRRADPROC_NO_INTERPOLATE_SUNRISE_SUNSET);
			irr.set_location(lat[locIndex], lon[locIndex], tz[locIndex]);
			
			irr.set_sky_model(1, 0.2);
			irr.set_beam_diffuse(dn, dh);
			irr.set_surface(0, tilt[locIndex], 180., 0., 0, 0.3);

			EXPECT_FALSE(irr.calc());

			// using HDKR model
			int sunup = 0;
			double solazi = 0, solzen = 0, solalt = 0;
			double ibeam, iskydiff, ignddiff;
			double aoi, stilt, sazi, rot, btd;

			irr.get_sun(&solazi, &solzen, &solalt, 0, 0, 0, &sunup, 0, 0, 0);
			irr.get_angles(&aoi, &stilt, &sazi, &rot, &btd);
			irr.get_poa(&ibeam, &iskydiff, &ignddiff, 0, 0, 0);


			pvinput_t in(ibeam, iskydiff, ignddiff, 0, 0,
				T_cell, 0.0, 0.0, 0.0, 0.0,
				solzen, aoi, elev[locIndex],
				stilt, sazi,
				((double)hour) + minute / 60.0,
				Irradiance_IO::DN_DF, false);
			pvoutput_t out;


			if (currentModule.name != name) {
				for (size_t n = 0; n < mmVector.size(); n++) {
					if (mmVector[n].name == name) currentModule = mmVector[n];
				}
			}
			cec6par_module_t solver;

			solver.a = currentModule.A;
			solver.Rs = currentModule.Rs;
			solver.Rsh = currentModule.Rsh;
			solver.Io = currentModule.Io;
			solver.Il = currentModule.Il;
			solver.Adj = currentModule.Adj;
			solver.alpha_isc = currentModule.alphaIsc;
			solver.beta_voc = currentModule.betaVoc;
			solver.Voc = currentModule.Voc;
			solver.Imp = currentModule.Imp;
			solver.Isc = currentModule.Isc;
			

			double voltage = -1;
			solver(in, T_cell, voltage, out);

			outputFile << fileLines[i+1] << "," << out.Power << ",";

			// using measured POA
			pvinput_t in2(0., 0., 0., 0, POA_meas,
				T_cell, 0.0, 0.0, 0.0, 0.0,
				solzen, aoi, elev[locIndex],
				stilt, sazi,
				((double)hour) + minute / 60.0,
				3, true);
			solver(in2, T_cell, -1, out);
			outputFile << out.Power << "\n";
		}
		outputFile.close();
	}

}
