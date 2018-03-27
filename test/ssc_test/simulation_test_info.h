#ifndef __TestInfo__
#define __TestInfo__

#include <string>

enum {
	STR,
	NUM,
	ARR,
	MAT
};

struct TestInfo{
	const char* sscVarName;		///< same name as SSC variable
	unsigned int dataType;		///< 0: string, 1: number, 2: array, 3: matrix

	const char* values;			///< comma-separated for array of values
	size_t length = 1;			///< length of array/matrix
	size_t width = 1;			///< width of matrix
};

// enum for test_types: equal, near (approx equal), greater than, less than, bool, cmod error
enum {
	EQ,
	NR,
	GT,
	LT,
	BL,
	ERR
};

struct TestResult {
	//TestInfo* test_values;

	const char* sscVarName;
	unsigned int testType;	
	double expectedResult;	///< expected test results per value
	double errorBound;			///< percent error allowed
};

static void modifyDefaults(TestInfo defaults, TestInfo specificCase) {
	//int x;
};

class SimulationTestTable {
public:
	SimulationTestTable(std::string computeModule, TestInfo* I, int nInfo, TestResult* R, int nRes) {
		computeModuleType = computeModule;
		info = I;
		result = R;
		nI = nInfo;
		nR = nRes;
	}
	int getNumInfo() { return nI; }
	int getNumResult() { return nR; }
	TestInfo* getInfo() { return info; }
	TestResult* getResult() { return result; }
private:
	std::string computeModuleType;
	TestInfo* info;
	TestResult* result;
	int nI, nR;
};

#endif
