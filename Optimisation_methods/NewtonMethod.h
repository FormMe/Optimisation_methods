#pragma once
#include "Solver.h"

class NewtonMethod: public Solver
{
public:
	NewtonMethod(ifstream &fin);
	~NewtonMethod(){};

	virtual Vertex Calc(func _f);

protected:
	Vertex x1;
	Vertex x2;
	Vertex x3;
	Vertex grad;
	vector<vector<double>> H;
	vector<vector<double>> H1;
	double lambda;

	virtual void Hessian();
	void Grad();
	void Inversion();
};
