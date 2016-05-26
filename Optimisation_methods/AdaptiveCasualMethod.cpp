#include "AdaptiveCasualMethod.h"
#include <random>


Vertex AdaptiveCasualMethod::Calc(func _f, const Vertex& _x)
{
	f = _f; x = _x;
	auto prevF = f(x.vec) + 1;
	funcCnt += 2;
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-1., 1.);
	auto alpha = _alpha;
	auto beta = _beta;
	auto step = _step;
	auto minStep = _minStep;
	auto quit = false;
	for (auto i = 0; i < M && !quit; i++)
	{
		for (auto j = 0; j < trials; ++j)
		{
			for (auto &coord : direct.vec)
				coord = dis(gen);
			auto fx = f(x.vec);
			auto y = x + (direct / direct.norm()) * step;
			if (f((y).vec) >= fx) continue;
			auto z = x + (y - x)*alpha;
			if (f(z.vec) >= fx) continue;
			x = z;
			step *= alpha;
			break;
		}
		quit = step <= minStep;
		if (!quit) step *= beta;
	}
	return x;
}
