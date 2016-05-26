#pragma once
#include "Solver.h"

class AdaptiveCasualMethod : public Solver
{
public:
	AdaptiveCasualMethod(ifstream &fin) :
		Solver(fin)
	{
		fin >> trials >> _step >> _minStep >> _alpha >> _beta;
		direct = Vertex(N);
	}

	virtual Vertex Calc(func _f, const Vertex &x) override;

private:
	int trials;
	double _step;
	double _minStep;
	double _alpha;
	double _beta;
	Vertex direct;


	bool QuitCase();
};

