#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "lib_iec61853_test.h"




TEST_F(IEC61853Test, ParameterEstimateWithIECSolverTest) {
	for (size_t m = 0; m < 20; m++) {
		//[IRR, TC, PMP, VMP, VOC, ISC]
		std::vector<double> testData(108);
		for (size_t n = 0; n < 18; n++) {
			testData[6 * n] = mmVector[m * 18 + n].irradiance;
			testData[6 * n + 1] = mmVector[m * 18 + n].temp;
			testData[6 * n + 2] = mmVector[m * 18 + n].Pm;
			testData[6 * n + 3] = mmVector[m * 18 + n].Vmp;
			testData[6 * n + 4] = mmVector[m * 18 + n].Voc;
			testData[6 * n + 5] = mmVector[m * 18 + n].Isc;
		}
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

void mychoice_switch(int choice) {
	switch (choice) {
	case -10:
		printf("You are being very negative!");
	case 0:
		printf("The invention of zero was brilliant!");
	default:
		printf("Meh!");
	}
}

TEST_F(IEC61215Test, unsolvedModules) {
	// Uncertainty for crystalline silicon modules: 
	// Pm = ± 2.8%,  Isc = ± 2.3%, Imp = ± 2.3%, Voc = ± 0.3%, Vmp = ± 0.7%
	std::ofstream file;
	file.open("C:/Users/dguittet/Documents/IEC 61853 Modeling/Data For Validating Models/loopcheck.csv");
	EXPECT_TRUE(file.is_open());
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

		int intervals = 10;
		double dPm = mmVector[i].Pm * .05 / (double)(intervals/2.);
		double dIsc = mmVector[i].Isc * 0.05 / (double)(intervals / 2.);
		double dImp = mmVector[i].Imp * 0.05 / (double)(intervals / 2.);
		double dVoc = mmVector[i].Voc * 0.05 / (double)(intervals/2);
		double dVmp = mmVector[i].Vmp * 0.05 / (double)(intervals/2);
		double Pm0 = mmVector[i].Pm;
		double Isc0 = mmVector[i].Isc;
		double Voc0 = mmVector[i].Voc;
		double Vmp0 = mmVector[i].Vmp;
		double Imp0 = mmVector[i].Imp;
		mmVector[i].Isc -= mmVector[i].Isc * 0.05;
		mmVector[i].Imp -= mmVector[i].Imp * 0.05;
		mmVector[i].Voc -= mmVector[i].Voc * 0.05;
		mmVector[i].Vmp -= mmVector[i].Vmp * 0.05;
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
						if (abs((mmVector[i].Pm - Pm0) / Pm0) > 0.05) continue;
						if (mmVector[i].Vmp > mmVector[i].Voc) continue;
						int tech_id = mmVector[i].type;
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
						double tol = 1e-10;
						while (attempts < 8 && mmVector[i].solved != 0) {
							//mmVector[i].solved = m.solve_with_sanity_and_heuristics<double>(300, tol);
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
						file << mmVector[i].Imp << ", " << mmVector[i].Vmp << ", ";
						file << mmVector[i].Isc << ", " << mmVector[i].Voc << "," << std::endl;
					
						if (mmVector[i].solved == 0) break;	
					}
					if (mmVector[i].solved == 0) break;
					mmVector[i].Vmp = Vmp0 - Vmp0 * 0.05;
					
				}
				if (mmVector[i].solved == 0) break;
				mmVector[i].Voc = Voc0 - Voc0 * 0.05;
			}
			mmVector[i].Imp = Imp0 - Imp0 * 0.05;
			if (mmVector[i].solved == 0) break;
		}
		if (mmVector[i].solved != 0) {
			mmVector[i].Io = Io_best;
			mmVector[i].Il = Il_best;
			mmVector[i].Rs = Rs_best;
			mmVector[i].Rsh = Rsh_best;
			mmVector[i].Adj = Adj_best;
			mmVector[i].A = a_best;
			mmVector[i].solved = solved_best;
		}
	}
	//mmVectorToCSV();
	file.close();
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
		int tech_id = mmVector[i].type;
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

TEST_F(IEC61215Test, testMeasurements) {
	
	getTestMeasurementData();

	std::string outputFileName = "C:\\Users\\dguittet\\Documents\\IEC 61853 Modeling\\Data For Validating Models\\singleDiodeResults.csv";
	std::ofstream outputFile;
	outputFile.open(outputFileName);

	moduleMeasurements currentModule;
	std::vector<double> powerPredicted;
	for (size_t i = 0; i < testMeasurements.size(); i++) {
		std::string name = testMeasurements[i][0].c_str();
		double Geff = atof(testMeasurements[i][2].c_str());
		double T_cell = atof(testMeasurements[i][4].c_str());

		double a, rs, rsh, io, il, adj;
		bool solved = false;
		if (currentModule.name != name) {
			for (size_t n = 0; n < mmVector.size(); n++) {
				if (mmVector[n].name == name) currentModule = mmVector[n];
			}
		}
		a = currentModule.A;
		rs = currentModule.Rs;
		rsh = currentModule.Rsh;
		io = currentModule.Io;
		il = currentModule.Il;
		adj = currentModule.Adj;

		if (!solved) continue;
		
		double eg0 = 1.12;
		double Tc_ref = (25 + 273.15);
		double KB = 8.618e-5;

		T_cell = T_cell + 273.15; // want cell temp in kelvin
		double muIsc = currentModule.alphaIsc * (1 - currentModule.Adj / 100);
		// calculation of IL and IO at operating conditions
		double IL_oper = Geff / 1000. * (il + muIsc * (T_cell - 298.15));
		if (IL_oper < 0.0) IL_oper = 0.0;

		double EG = eg0 * (1 - 0.0002677*(T_cell - (25 + 273.15)));
		double IO_oper = currentModule.Io * pow(T_cell / Tc_ref, 3) * exp(1 / KB * (eg0 / Tc_ref - EG / T_cell));
		double A_oper = a * T_cell / Tc_ref;
		double Rsh_oper = currentModule.Rsh * (1000. / Geff);

		double V_oc = openvoltage_5par(currentModule.Voc, A_oper, IL_oper, IO_oper, Rsh_oper);
		double I_sc = IL_oper / (1 + currentModule.Rs / Rsh_oper);


		double V, I;
		maxpower_5par(100, a, il, io, rs, rsh, &V, &I);
		powerPredicted.push_back(V*I);
		outputFile << V * I << "\n";
	}
	outputFile.close();

}

