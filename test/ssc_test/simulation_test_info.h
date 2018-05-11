#ifndef __TestInfo__
#define __TestInfo__

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
	TF,
	ERR
};

struct TestResult {
	//TestInfo* test_values;

	const char* sscVarName;
	unsigned int testType;	
	double expectedResult;		///< expected test results per value
	double errorBound;			///< percent error allowed
};

class SimulationTestTable {
public:
	const char* name;
	SimulationTestTable(const char* cmodType, const char* testName, TestInfo* I, int nInfo, TestResult* R, int nRes) {
		name = testName;
		computeModuleType = cmodType;
		info = I;
		result = R;
		nI = nInfo;
		nR = nRes;
	}
	SimulationTestTable(const SimulationTestTable& first) {
		name = first.name;
		computeModuleType = first.computeModuleType;
		info = first.info;
		result = first.result;
		nI = first.nI;
		nR = first.nR;
	}
	const char* getCMODType() { return computeModuleType; }
	int getNumInfo() { return nI; }
	int getNumResult() { return nR; }
	TestInfo* getInfo() { return info; }
	TestResult* getResult() { return result; }
	
	void swap(SimulationTestTable& other) {
		using std::swap;
		SimulationTestTable first = *this;
		swap(first.name, other.name);
		swap(first.computeModuleType, other.computeModuleType);
		swap(first.info, other.info);
		swap(first.result, other.result);
		swap(first.nI, other.nI);
		swap(first.nR, other.nR);
	}
	SimulationTestTable& operator=(SimulationTestTable other) {
		swap(other);
		return *this;
	}

protected:
	const char* computeModuleType;
	TestInfo* info;
	TestResult* result;
	int nI, nR;
};

void modifyTestInfo(SimulationTestTable* defaults, SimulationTestTable* specificCase) {
	int x;
};

#endif
