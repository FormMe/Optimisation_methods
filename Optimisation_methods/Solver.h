#pragma once
#include "Vertex.h"
#include "LU_eigenvalues.h"

class Solver
{
public:
	Solver(): eps(1e-12), h(1e-10), h1(1e-10), M(10000) {};
	Solver(string filename);
	virtual ~Solver() {};
	virtual Vertex Calc(func _f, Vertex &_x) { return{}; };

protected:
	LU_eigenvalues _eigenvalues;
	func f;
	Vertex x;
	Vertex S;
	double eps;
	double h;
	double h1;
	int N;
	int M;

	double GSS(Vertex &_x, Vertex &_S);

private:
	pair<double, double> FindInterval(double x0, double d, Vertex &_x, Vertex &_S);
};

