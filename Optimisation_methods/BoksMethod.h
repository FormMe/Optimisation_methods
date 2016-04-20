#pragma once
#include "Solver.h"
#include <random>

class BoksMethod : Solver
{
public:
	BoksMethod(ifstream &fin);
	~BoksMethod() {};

	Vertex Calc(func _f);

private:
	int K;
	vector<Vertex> cmplx;
	vector<double> F;

	int max_ind;
	double alpha;
	
	Vertex L;
	Vertex R;
	Vertex average;
	
	void InitCmplx();
	void Average();
	bool QuitCase();
	void CorrectVertex(Vertex &v);
};

