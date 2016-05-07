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
	Vertex x1;
	Vertex grad;
	Vertex grad1;
	Vertex S;
	double lambda;
	double w;

	void Grad();
	void PolakRibiere();
	void FletcherReeves();



};

