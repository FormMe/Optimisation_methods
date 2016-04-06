#pragma once
#include "Vertex.h"
#include "Solver.h"

class NonlinearConjugateGradientMethod : Solver
{
public:
	NonlinearConjugateGradientMethod(string filename) : Solver(filename) {};
	NonlinearConjugateGradientMethod() : Solver() {};
	~NonlinearConjugateGradientMethod();

	Vertex NCGM(func _f, Vertex _x);

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

