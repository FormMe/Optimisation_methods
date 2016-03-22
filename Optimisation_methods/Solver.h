#pragma once
#include "DownhillSimplexMethod.h"
#include "NonlinearConjugateGradientMethod.h"


class Solver : public NonlinearConjugateGradientMethod
{
public:
	Solver();
	~Solver();
};

