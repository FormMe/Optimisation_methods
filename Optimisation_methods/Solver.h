#pragma once
#include "Vertex.h"
#include "LU_eigenvalues.h"

class Solver
{
public:
	Solver(ifstream &fin);
	Solver(
		double eps,
		double eps1,
		double l,
		double lStep,
		double h,
		double h1,
		int N,
		int M,
		int M1);

	virtual Vertex Calc(func _f, const Vertex &x) = 0;


	int GetFuncCnt();

protected:
	LU_eigenvalues _eigenvalues;
	func f;
	Vertex x;
	Vertex x1;
	Vertex grad;
	Vertex prevX;
	Vertex S;
	Vertex L;
	Vertex R;
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

	double Fibbonachi(Vertex &_x, Vertex &_S);
	double GSS(Vertex &_x, Vertex &_S);
	void CorrectVertex(Vertex &v);
	void Grad();

private:
	pair<double, double> FindInterval(double x0, double d, Vertex &_x, Vertex &_S);
};

