#pragma once
#include "Vertex.h"
#include "GoldenSectionSearch.h"

class NonlinearConjugateGradientMethod
{
public:
	NonlinearConjugateGradientMethod();
	~NonlinearConjugateGradientMethod();

	Vertex NCGM(func _f, Vertex _x);

private:
	GoldenSectionSearch _gss;
	func f;
	Vertex x;
	Vertex x1;
	Vertex grad;
	Vertex grad1;
	Vertex S;
	double lambda;
	double w;
	double eps;
	double h;
	int N;

	void Grad();
	void PolakRibiere();
	void FletcherReeves();



};

