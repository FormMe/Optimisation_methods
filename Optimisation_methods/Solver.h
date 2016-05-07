#pragma once
#include "Vertex.h"
#include "LU_eigenvalues.h"

class Solver
{
public:
	Solver(ifstream &fin);
	virtual Vertex Calc(func _f, const Vertex &x) { return{}; };

	int GetFuncCnt();

protected:
	LU_eigenvalues _eigenvalues;
	func f;
	Vertex x;
	Vertex prevX;
	Vertex S;
	double eps;
	double eps1;
	double l;
	double lStep;
	double h;
	double h1;
	int N;
	int M;
	int M1;
	int funcCnt;

	double GSS(Vertex &_x, Vertex &_S);

private:
	pair<double, double> FindInterval(double x0, double d, Vertex &_x, Vertex &_S);
};

