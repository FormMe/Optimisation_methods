#pragma once
#include "Solver.h"
#include <random>

class StaticalGradMethod: public Solver
{
public:
	StaticalGradMethod(ifstream &fin) : 
		Solver(fin),
		dF(Vertex(N))
	{
		fin >> trials >> step;
		directions = vector<Vertex>(trials, Vertex(N));
	}

	virtual Vertex Calc(func _f, const Vertex &x);

private:
	int trials;
	double step;
	double lambda;
	vector<Vertex> directions;
	Vertex dF;

	void GenerateDirections();
	void Calculate_dF();
};

