#include "PenaltyMethod.h"
#include <numeric>

PenaltyMethod::PenaltyMethod(ifstream &fin, const func &_f, const vector<func> &_g, Solver *_s) :
	s(_s),
	f(_f),
	g(_g)
{
	fin >> N >> M >> r >> C >> penalty_eps;
	x = Vertex(N);
	for (auto &x0 : x.vec)
		fin >> x0;
}

PenaltyMethod::PenaltyMethod(double _C, double _r, double _penalty_eps, int _M, Vertex _x, const func& _f, const vector<func>& _g, Solver* _s) :
	C(_C),
	r(_r),
	penalty_eps(_penalty_eps),
	M(_M),
	x(_x),
	N(x.vec.size()),
	s(_s),
	f(_f),
	g(_g) {	}



Vertex PenaltyMethod::Calculate()
{

	auto quit = false;
	for (auto k = 0; k < M + 100 && !quit; k++)
	{
		////רענאפû

		auto P = [&](vector<double> vec)
		{
			return  r / 2 * accumulate(g.begin(), g.end(), 0.0,
				[&vec](double a, const func &b)
			{
				auto fb = b(vec);
				return a + max(0.0, fb) * fb;
			});
		};

		////באנüונû

		//auto P = [&](vector<double> vec)
		//{
		//	auto sum = accumulate(g.begin(), g.end(), 0.0,
		//		[&vec](double a, const func &b)
		//	{
		//		if (a == DBL_MAX) return a;
		//		auto lim = b(vec);
		//		return lim < 0 ? (a + 1 / lim) : DBL_MAX;
		//	});
		//	return sum == DBL_MAX ? sum : -r * sum;
		//};

		//auto P = [&](vector<double> vec)
		//{
		//	auto sum = accumulate(g.begin(), g.end(), 0.0,
		//		[&vec](double a, const func &b)
		//	{
		//		if (a == DBL_MAX) return a;
		//		auto lim = b(vec);
		//		return lim < 0 ? (a + log(-lim)) : DBL_MAX;
		//	});
		//	return sum == DBL_MAX ? sum : -r * sum;
		//};

		auto F = [&](const vector<double> &vec)
		{
			auto p = P(vec);
			return  p == DBL_MAX ? p : f(vec) + p;
		};
		x = s->Calc(F, x);
		auto px = P(x.vec);
		quit = abs(px) <= penalty_eps;
		r *= C;
	}
	return x;
}

int PenaltyMethod::GetFuncCount()
{
	return  s->GetFuncCnt();
}

Vertex PenaltyMethod::Calc()
{
	return Calculate();
}

Vertex PenaltyMethod::Calc(Vertex _x)
{
	x = _x;
	return Calculate();
}
