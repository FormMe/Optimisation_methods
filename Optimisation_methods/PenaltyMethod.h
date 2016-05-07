#pragma once
#include "Vertex.h"
#include "Solver.h"
#include <numeric>

class PenaltyMethod
{
private:
	Solver *s;
	Vertex x;
	func f;
	vector<func> g;

	double C;
	double r;
	double penalty_eps;
	int M;
	int N;

public:
	PenaltyMethod(ifstream &fin, const func &_f, const vector<func> &_g, Solver *_s);

	Vertex Calc();
};

