#pragma once
#include "Solver.h"
#include <random>

class BoksMethod : Solver
{
public:
	BoksMethod(string filename) : Solver(filename) {};
	BoksMethod() : Solver() {};
	~BoksMethod() {};

	Vertex Calc(func _f, Vertex &v, Vertex &L, Vertex &R, double alpha);

private:
	int K;
	vector<Vertex> cmplx;
	vector<double> F;

	int max_ind;
	
	Vertex L;
	Vertex R;
	Vertex average;
	
	void Initialization(Vertex &v);
	void Average();
	bool QuitCase();
	void CorrectVertex(Vertex &v);
};

