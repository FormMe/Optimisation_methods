#pragma once
#include "Solver.h"

class RosenbrockMethod : Solver
{
public:
	RosenbrockMethod(string filename) : Solver(filename) {};
	RosenbrockMethod() : Solver(){};
	~RosenbrockMethod(){};
	virtual Vertex Calc(func _f, Vertex &_x);

private:
	Vertex x1;
	vector<Vertex> S;
	vector<Vertex> A;
	vector<double> lambda;

	void MinDirections();
	void FindDirectios();
	void GramSchmidtProcess();
	void PalmerProcess();

};

