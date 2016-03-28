#pragma once
#include "DownhillSimplexMethod.h"
#include "NonlinearConjugateGradientMethod.h"
#include "NewtonMethod.h"

class Solver : public NonlinearConjugateGradientMethod,
			public NewtonMethod
{
public:
	Solver();
	~Solver();
};

