#ifndef __lib_iec61853_test__
#define __lib_iec61853_test__

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <lib_util.h>
#include <6par_solve.h>
#include <lib_iec61853.h>
#include <lib_cec6par.h>

#include <gtest/gtest.h>

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
	double Il = 0.0;
	double Io = 0.0;
	double Rs = 0.0;
	double Rsh = 0.0;
	double A = 0.0;
	double Adj = 0.0;
	double alphaIsc = 0.0;
	double betaVoc = 0.0;
	double gammaPmp = 0.0;
	int solved = -1;

	std::vector<std::vector<double>> currentAt10Voltages; // divided Voc by 10 for step size
};

class IEC61853UpdateSolver {
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
};

class IEC61853Test : public ::testing::Test {
protected:
	std::string moduleMeasurementsFile;
	IEC61853UpdateSolver solver;
	std::vector<moduleMeasurements> mmVector;
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
				file << mmVector[r].name << "," << mmVector[r].Voc << "," << mmVector[r].temp << "," << mmVector[r].irradiance << ", ";
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
			"Il_iec",  "Io_iec", "Rs_iec", "Rsh_iec", "a_iec", "Adj_cec",
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
			file << mmVector[r].name << "," << mmVector[r].temp << "," << mmVector[r].irradiance << ", ";
			file << mmVector[r].Isc << "," << mmVector[r].Voc << "," << mmVector[r].Imp << "," << mmVector[r].Vmp << "," << mmVector[r].Pm << "," << mmVector[r].nSer << ", ";
			file << mmVector[r].Il << ", " << mmVector[r].Io << ", " << mmVector[r].Rs << ", " << mmVector[r].Rsh << ", " << mmVector[r].A << ", " << mmVector[r].Adj << ",";
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
			size_t start = pos + 1;
			std::vector<double> v;
			for (size_t i = 0; i < 14; i++) {
				size_t pos = str.find(delimiter, start);
				std::string token = str.substr(start, pos - start);
				double val = atof(token.c_str());
				v.push_back(val);
				start = pos + 1;
			}
			mm.temp = v[0];
			mm.irradiance = v[1];
			mm.Isc = v[2];
			mm.Voc = v[3];
			mm.Imp = v[4];
			mm.Vmp = v[5];
			mm.Pm = v[6];
			mm.nSer = (int)v[7];
			mm.Il = v[8];
			mm.Io = v[9];
			mm.Rs = v[10];
			mm.Rsh = v[11];
			mm.A = v[12];
			mm.solved = v[13];

			mmVector.push_back(mm);
		}
		file.close();
	}
};

class IEC61215Test : public ::testing::Test {
protected:
	std::string moduleMeasurementsFile;
	std::vector<moduleMeasurements> mmVector;

	std::string testMeasurementsFiles;
	std::vector<std::vector<std::string>> testMeasurements;
	double e = 0.01;
public:
	void mmVectorToCSV() {
		std::ofstream file;
		file.open(moduleMeasurementsFile);
		EXPECT_TRUE(file.is_open());

		std::vector<std::string> colHeaders = { "Name","alphaIsc",	"betaVoc", "gammaPmp", "Nser", "Temp", "Irradiance",
			"Isc", "Voc", "Imp", "Vmp", "Pm",
			"Il",  "Io", "Rs", "Rsh", "a", "Adj",
			"Solved" };
		for (size_t i = 0; i < colHeaders.size(); i++) {
			file << colHeaders[i] << ",";
		}
		file << "\n";
		for (size_t r = 0; r < mmVector.size(); r++) {
			file << mmVector[r].name << "," << mmVector[r].alphaIsc << "," << mmVector[r].betaVoc << "," << mmVector[r].gammaPmp << "," << mmVector[r].nSer << ", ";
			file << mmVector[r].temp << "," << mmVector[r].irradiance << ", ";
			file << mmVector[r].Isc << "," << mmVector[r].Voc << "," << mmVector[r].Imp << "," << mmVector[r].Vmp << "," << mmVector[r].Pm << "," ;
			file << mmVector[r].Il << ", " << mmVector[r].Io << ", " << mmVector[r].Rs << ", " << mmVector[r].Rsh << ", " << mmVector[r].A << ", " << mmVector[r].Adj << ",";
			file << mmVector[r].solved << ", \n";
		}
		file.close();
	}

	void getTestMeasurementData() {
		testMeasurementsFiles = "C:\\Users\\dguittet\\Documents\\IEC 61853 Modeling\\Data For Validating Models\\testMeasurements.csv";
		std::ifstream file;
		file.open(testMeasurementsFiles);
		EXPECT_TRUE(file.is_open());
		std::string str;
		std::string delimiter = ",";
		getline(file, str); // header row
		std::vector<double> powerPredicted;
		while (getline(file, str)) {
			size_t pos = str.find(delimiter, 0);
			std::string token = str.substr(0, pos);
			size_t start = pos + 1;
			std::vector<std::string> v;
			v.push_back(token);
			for (size_t i = 0; i < 8; i++) {
				size_t pos = str.find(delimiter, start);
				std::string token = str.substr(start, pos - start);
				v.push_back(token);
				start = pos + 1;
			}
			testMeasurements.push_back(v);
		}
		file.close();
	}

	void SetUp() {
		moduleMeasurementsFile = "C:/Users/dguittet/Documents/IEC 61853 Modeling/Data For Validating Models/moduleMeasurementsLDRD.csv";
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
			size_t start = pos + 1;
			std::vector<double> v;
			for (size_t i = 0; i < 18; i++) {
				size_t pos = str.find(delimiter, start);
				std::string token = str.substr(start, pos - start);
				double val = atof(token.c_str());
				v.push_back(val);
				start = pos + 1;
			}
			mm.temp = 25.;
			mm.irradiance = 1000.;
			mm.alphaIsc = v[0];
			mm.betaVoc = v[1];
			mm.gammaPmp = v[2];
			mm.nSer = (int)v[3];
			mm.Isc = v[6];
			mm.Voc = v[7];
			mm.Imp = v[8];
			mm.Vmp = v[9];
			mm.Pm = v[10];
			mm.Il = v[11];
			mm.Io = v[12];
			mm.Rs = v[13];
			mm.Rsh = v[14];
			mm.A = v[15];
			mm.Adj = v[16];
			mm.solved = v[17];


			mmVector.push_back(mm);
		}
		file.close();
	}
};
#endif