#ifndef _LIB_SINGLETON_TEST_H_
#define _LIB_SINGLETON_TEST_H_

#include <gtest/gtest.h>

#include "lib_singleton.h"

class TestClass
{
public:
	
	void SetData(int x_in) {
		x = x_in;
	}
	int GetData() {
		return x;
	}


private:
	TestClass() {};
	~TestClass() {};
	
	int x;
};

typedef SingletonHolder<TestClass, CreateUsingNew> SingleTestClass;

class SingletonTest : public ::testing::Test
{
public:

	void SetUp()
	{
		T = SingleTestClass::Instance();
		T.SetData(10);
	}
	void TearDown() {};
	 
protected:
	TestClass T;
};

TEST_F(SingletonTest, Create)
{
	EXPECT_EQ(T.GetData(), 10);
}

#endif
