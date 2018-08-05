#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <lib_util.h>
#include <6par_solve.h>
#include <lib_iec61853.h>

struct moduleMeasurements {
	std::string name;
	int type = -1;
	double temp;
	double irradiance;
	int nSer;
	double Isc;
	double Voc;
	double Imp;
	double Vmp;
	double Pm;
	double Il_iec = 0.0;
	double Io_iec = 0.0;
	double Rs_iec = 0.0;
	double Rsh_iec = 0.0;
	double A_iec = 0.0;
	double Il_cec = -1;
	double Io_cec = -1;
	double Rs_cec = -1;
	double Rsh_cec = -1;
	double A_cec = -1;
	double Adj_cec = -1;
	int solved = -1;

	std::vector<std::vector<double>> currentAt10Voltages; // divided Voc by 10 for step size
};

class IEC61853UpdateSolver{
public:
	double Il, Io, Rs, Rsh, A;	// values of current iteration
	double minIl, minIo, minRs, minRsh, minA;
	int nSer;
	double uncertainty;			// measured values uncertainty
	double tolerance;			// error tolerance for single diode model equations

	iec61853_module_t modIEC;
	module6par mod6p;

	void guessIEC(moduleMeasurements mm) {
		double Rs_scale, Rsh_scale;
		switch (mm.type)
		{
		case iec61853_module_t::Amorphous:
			Rs_scale = 0.59;
			Rsh_scale = 0.922;
			break;
		case iec61853_module_t::CdTe:
			Rs_scale = 0.46;
			Rsh_scale = 1.11;
			break;
		case iec61853_module_t::CIGS:
			Rs_scale = 0.55;
			Rsh_scale = 1.22;
			break;
		case iec61853_module_t::CIS:
			Rs_scale = 0.61;
			Rsh_scale = 1.07;
			break;
		case iec61853_module_t::monoSi:
			Rs_scale = 0.32;
			Rsh_scale = 4.92;
			break;
		case iec61853_module_t::multiSi:
		default:
			Rs_scale = 0.34;
			Rsh_scale = 5.36;
			break;
		}

		double Rs_ref0 = Rs_scale * (mm.Voc - mm.Vmp) / mm.Imp;

		if (Rs_ref0 < 0.02) Rs_ref0 = 0.02;
		if (Rs_ref0 > 60) Rs_ref0 = 60;

		double Rsh_ref0 = Rsh_scale * mm.Voc / (mm.Isc - mm.Imp);

		A = 3;
		Il = 0.95*mm.Isc;
		Rsh = Rsh_ref0;
		Io = (Il - mm.Voc / Rsh) / (exp(mm.Voc / A) - 1);
		Rs = Rs_ref0;

	}

	void setInitialValues(double Il0, double Io0, double Rs0, double Rsh0, double a0) {

	}

	void updateMinErrorSol(double error) {

	}

	double calculateError() {

	}
};


class IEC61853Test : public ::testing::Test {
protected:
	std::string moduleMeasurementsFile;
	size_t nColsmoduleMeasurementsFile = 16;
	std::vector<moduleMeasurements> mmVector;
	IEC61853UpdateSolver solver;
	double e = 0.01;
public:
	void mmCurrentToCSV() {
		std::ofstream file;
		file.open("C:/Users/dguittet/Documents/IEC 61853 Modeling/Data For Validating Models/moduleIVCurves.csv");
		EXPECT_TRUE(file.is_open());

		std::vector<std::string> colHeaders = { "Name", "Voc", "Temp", "Irradiance", "V", "I" };
		for (size_t i = 0; i < colHeaders.size(); i++) {
			file << colHeaders[i] << ",";
		}
		file << "\n";

		for (size_t r = 0; r < mmVector.size(); r++) {
			for (size_t i = 0; i < 10; i++) {
				file << mmVector[r].name << "," << mmVector[r].Voc<< "," << mmVector[r].temp << "," << mmVector[r].irradiance << ", ";
				file << mmVector[r].currentAt10Voltages[i][0] << ", " << mmVector[r].currentAt10Voltages[i][1] << ", ";
				file << "\n";
			}

		}
		file.close();
	}
	void mmVectorToCSV() {
		std::ofstream file;
		file.open(moduleMeasurementsFile);
		EXPECT_TRUE(file.is_open());

		std::vector<std::string> colHeaders = { "Name", "Temp", "Irradiance", "Isc", "Voc", "Imp", "Vmp", "Pm", "Nser",
			"Il_iec",  "Io_iec", "Rs_iec", "Rsh_iec", "a_iec",
			"Sanity err" };
			//"Solved", "Il_cec", "Io_cec", "Rs_cec", "Rsh_cec", "a_cec", "Adj_cec"};
		for (size_t i = 0; i < colHeaders.size(); i++) {
			file << colHeaders[i] << ",";
		}
		file << "V step" << ",";
		for (size_t i = 0; i < 10; i++) {
			file << "I_" << i << ",";
		}
		file << "\n";


		for (size_t r = 0; r < mmVector.size(); r++) {
			file << mmVector[r].name << "," << mmVector[r].temp << "," << mmVector[r].irradiance << ", " ;
			file << mmVector[r].Isc << "," << mmVector[r].Voc << "," << mmVector[r].Imp << "," << mmVector[r].Vmp << "," << mmVector[r].Pm << "," << mmVector[r].nSer << ", ";
			file << mmVector[r].Il_iec << ", " << mmVector[r].Io_iec << ", " << mmVector[r].Rs_iec << ", " << mmVector[r].Rsh_iec << ", " << mmVector[r].A_iec << ", ";
			file << mmVector[r].solved << ", ";
			/*if (mmVector[r].solved == 0) {
				file << mmVector[r].Il_cec << ", " << mmVector[r].Io_cec << ", " << mmVector[r].Rs_cec << ", " << mmVector[r].Rsh_cec << ", " << mmVector[r].A_cec << ", " << mmVector[r].Adj_cec << ",";
			}*/
			file << mmVector[r].currentAt10Voltages[1][0] << ", ";
			for (size_t i = 0; i < 10; i++) {
				file << mmVector[r].currentAt10Voltages[i][1] << ", ";
			}
			file << "\n";

		}
		file.close();
	}
	void SetUp() {
		moduleMeasurementsFile = "C:/Users/dguittet/Documents/IEC 61853 Modeling/Data For Validating Models/moduleMeasurements.csv";
		std::ifstream file;
		file.open(moduleMeasurementsFile);
		EXPECT_TRUE(file.is_open());
		std::string str;
		std::string delimiter = ",";
		getline(file, str); // header row
		while (getline(file, str)) {
			moduleMeasurements mm;

			size_t pos = str.find(delimiter);
			mm.name = str.substr(0, pos);
			std::string typeName = mm.name.substr(0, 3);

			if (typeName == "xSi") {
				mm.type = module6par::monoSi;
			}
			else if (typeName == "mSi") {
				mm.type = module6par::multiSi;
			}
			else if (typeName == "CdT") {
				mm.type = module6par::CdTe;
			}
			else if (typeName == "CIG") {
				mm.type = module6par::CIGS;
			}
			else if (typeName == "HIT") {
				mm.type = module6par::Amorphous;
			}
			else if (typeName == "aSi") {
				mm.type = module6par::Amorphous;
			}
			else {
				EXPECT_TRUE(mm.type != -1); //error
			}
			size_t start = pos+1;
			std::vector<double> v;
			for (size_t i = 0; i < 8; i++) {
				size_t pos = str.find(delimiter, start);
				std::string token = str.substr(start, pos-start);
				double val = atof(token.c_str());
				v.push_back(val);
				start = pos+1;
			}
			mm.temp = v[0];
			mm.irradiance = v[1];
			mm.Isc = v[2];
			mm.Voc = v[3];
			mm.Imp = v[4];
			mm.Vmp = v[5];
			mm.Pm = v[6];
			mm.nSer = (int)v[7];
			mmVector.push_back(mm);
		}
		file.close();
	}
};


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
			mmVector[m * 18 + n].Il_iec = par.at(n, 0);
			mmVector[m * 18 + n].Io_iec = par.at(n, 1);
			mmVector[m * 18 + n].Rs_iec = par.at(n, 2);
			mmVector[m * 18 + n].Rsh_iec = par.at(n, 3);
			mmVector[m * 18 + n].A_iec = par.at(n, 4);
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
	
		m.Io = mm->Io_iec;
		m.Il = mm->Il_iec;
		m.Rs = mm->Rs_iec;
		m.Rsh = mm->Rsh_iec;
		m.a = mm->A_iec;
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

TEST_F(IEC61853Test, DISABLED_ParameterEstimationTest_lib_iec61853) {
	for (size_t i = 0; i < mmVector.size(); i++) {
		int tech_id = mmVector[i].type;
		double Vmp = mmVector[i].Vmp;
		double Imp = mmVector[i].Imp;
		double Voc = mmVector[i].Voc;
		double Isc = mmVector[i].Isc;
		double alpha = 0;
		double beta = 0;
		double gamma = 0;
		int nser = mmVector[i].nSer;

		// monoSi, multiSi, CdTe, CIS, CIGS, Amorphous 
		module6par m(tech_id, Vmp, Imp, Voc, Isc, beta, alpha, gamma, nser, 298.15);
		m.setUncertaintyPmp(5);
		//std::cout << "module " << i << ": ";
		mmVector[i].solved = m.solve_with_sanity_and_heuristics<double>(300, 1e-4);
		//std::cout << mmVector[i].solved << "\n";
	}

	mmVectorToCSV();

}