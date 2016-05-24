#include "StaticalGradMethod.h"

Vertex StaticalGradMethod::Calc(func _f, const Vertex& _x)
{
	f = _f; x = _x;
	auto prevF = f(x.vec) + 1;
	funcCnt+=2;
	for (auto i = 0; i < M && abs(prevF - f(x.vec)) > eps; i++)
	{
		GenerateDirections();
		Calculate_dF();
		auto S = dF / (-dF.norm());
		lambda = GSS(x, S);
		prevF = f(x.vec);
		funcCnt+=2;
		x = x + S*lambda;
	}
	return x;
}

void StaticalGradMethod::GenerateDirections()
{
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> dis(0, 100);
	for (auto &direct : directions)
		for (auto &coord : direct.vec)
			coord = dis(gen) % 2;
}

void StaticalGradMethod::Calculate_dF()
{
	dF = Vertex(N);
	auto fx = f(x.vec);
	funcCnt++;
	for (auto direct : directions)
	{
		funcCnt++;
		auto dfx = f((x + direct*step).vec) - fx;
		dF = dF + direct*dfx;
	}
}
