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
		minX(Vertex(N)),
		funcCnt(0)
	{
		
	};

	Vertex Calc(PenaltyMethod *solver, func f, int M);
	int GetFuncCount();

private:
	Vertex L;
	Vertex R;

	vector<pair<Vertex, Vertex>> sectors;

	int N;

	Vertex x;
	Vertex minX;


	int funcCnt;
	void GetRandomX();
};

