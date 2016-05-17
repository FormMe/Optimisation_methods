#pragma once
#include "PenaltyMethod.h"
class GlobalMinimizator
{
public:
	GlobalMinimizator(const Vertex &l, const Vertex &r) :
		L(l),
		R(r),
		N(R.vec.size()),
		x(Vertex(N)),
		minX(Vertex(N)) { };

	Vertex Calc(PenaltyMethod *solver, func f, int M);
private:
	Vertex L;
	Vertex R;
	Vertex x;
	Vertex minX;
	int N;

	void GetRandomX();
};

