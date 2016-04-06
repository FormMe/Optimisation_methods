#pragma once
#include "Solver.h"

class RosenbrockMethod : Solver
{
public:
	RosenbrockMethod(string filename) : Solver(filename) {};
	RosenbrockMethod() : Solver(){};
	~RosenbrockMethod();
	Vertex RM(func _f, Vertex &v);

private:
	Vertex x1;
	vector<Vertex> S;
	vector<Vertex> A;
	vector<double> lambda;

	void FindDirectios();
	void GramSchmidtProcess();
	void PalmerProcess();

};

