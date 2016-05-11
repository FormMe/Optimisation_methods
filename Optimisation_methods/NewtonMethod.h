#pragma once
#include "Solver.h"

class NewtonMethod: public Solver
{
public:
	NewtonMethod(ifstream &fin);
	~NewtonMethod(){};

	virtual Vertex Calc(func _f, const Vertex &_x);

protected:
	Vertex x2;
	Vertex x3;
	vector<vector<double>> H;
	vector<vector<double>> H1;
	double lambda;

	void Hessian();
	void Inversion();
};
