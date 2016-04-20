#pragma once
#include "Solver.h"

class NewtonMethod: public Solver
{
public:
	NewtonMethod(ifstream &fin);
	~NewtonMethod(){};

	virtual Vertex Calc(func _f);

private:
	Vertex x1;
	Vertex x2;
	Vertex x3;
	Vertex grad;
	vector<vector<double>> H;
	vector<vector<double>> H1;
	double lambda;

	void Grad();
	void Hessian();
	void Inversion();
	void Init();
};
