#pragma once
#include "DownhillSimplexMethod.h"
#include "NonlinearConjugateGradientMethod.h"


class Solver : public DownhillSimplexMethod,
	public NonlinearConjugateGradientMethod
{
public:
	Solver();
	~Solver();
};

