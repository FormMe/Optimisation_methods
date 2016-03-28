#pragma once
#include "Vertex.h"
#include "GoldenSectionSearch.h"
#include "LU_eigenvalues.h"

class NewtonMethod
{
public:
	NewtonMethod();
	~NewtonMethod();

	Vertex NM(func _f, Vertex &_x);

private:
	LU_eigenvalues _eigenvalues;
	GoldenSectionSearch _gss;
	func f;
	Vertex x;
	Vertex x1;
	Vertex x2;
	Vertex x3;
	Vertex grad;
	Vertex d;
	vector<vector<double>> H;
	vector<vector<double>> H1;
	double lambda;
	double eps;
	double h;
	double h1;
	int N;

	void Grad();
	void Hessian();
	void Inversion();
	void Init();
};
