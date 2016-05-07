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


Vertex PenaltyMethod::Calc()
{
	auto quit = false;
	for (auto k = 0; k < M + 100 && !quit; k++, r /= C)
	{
		////רענאפû

		//auto P = [&](vector<double> vec)
		//{
		//	return  r / 2 * accumulate(g.begin(), g.end(), 0.0,
		//		[&vec](double a, const func &b)
		//	{
		//		auto fb = b(vec);
		//		return a + max(0.0, fb) * fb;
		//	});
		//};

		////באנüונû

		auto P = [&](vector<double> vec)
		{
			return  -r * accumulate(g.begin(), g.end(), 0.0,
				[&vec](double a, const func &b)
			{
				auto lim = b(vec);
				return a + lim <= 0 ? 1 / lim : DBL_MAX;
			});
		};

		//auto P = [&](vector<double> vec)
		//{
		//	return  -r * accumulate(g.begin(), g.end(), 0.0,
		//		[&vec](double a, const func &b)
		//		{
		//			return a + log(-b(vec));
		//		});
		//};

		auto F = [&](const vector<double> &vec)
		{
			auto p = P(vec);
			return f(vec) + p;
		};

		x = s->Calc(F, x);
		auto px = P(x.vec);
		quit = abs(px) <= penalty_eps;
	}
	cout << s->GetFuncCnt() << endl;
	return x;
}