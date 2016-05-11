#pragma once
#include "Vertex.h"
#include "Solver.h"

class NonlinearConjugateGradientMethod : public Solver
{
public:
	NonlinearConjugateGradientMethod(ifstream &fin) ;
	~NonlinearConjugateGradientMethod(){};

	virtual Vertex Calc(func _f, const Vertex &_x);

private:
	Vertex grad1;
	Vertex S;
	double lambda;
	double w;

	void PolakRibiere();
	void FletcherReeves();



};

