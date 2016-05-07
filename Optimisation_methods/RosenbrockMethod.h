#pragma once
#include "Solver.h"

class RosenbrockMethod : public Solver
{
public:
	RosenbrockMethod(ifstream &fin);
	Vertex Calc(func _f, const Vertex &_x);

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

