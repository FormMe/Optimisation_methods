#pragma once
#include "NewtonMethod.h"

class LevenbergMarquardtMethod : public NewtonMethod
{
public:
	LevenbergMarquardtMethod(ifstream &fin);
	~LevenbergMarquardtMethod();

	virtual Vertex Calc(func _f);

private:
	void WeightH();
};

