#include "Solver.h"

Solver::Solver(ifstream &fin) : funcCnt(0)
{
	fin >> eps >> eps1 >> l >> lStep >> h >> h1 >> M >> M1 >> N;
	_eigenvalues = LU_eigenvalues(min(eps,eps1));
	x = Vertex(N);
	prevX = Vertex(N);
	S = Vertex(N);
	for (auto i = 0; i < N; ++i)
		fin >> x.vec[i];
}


pair<double, double> Solver::FindInterval(double lambda0, double d, Vertex &_x, Vertex &_S)
{
	double lambda1, lambda2;
	if (f((_x + _S*lambda0).vec) < f((_x + _S*(lambda0 + d)).vec))
		d = -d;
	do
	{
		lambda1 = lambda0 + d;
		d *= 2;
		lambda2 = lambda1 + d;
		lambda0 = lambda1;

	} while (f((_x + _S*lambda1).vec) > f((_x + _S*lambda2).vec));
	lambda0 -= d / 2;
	return make_pair(min(lambda0, lambda2), max(lambda0, lambda2));
}


double Solver::GSS(Vertex &_x, Vertex &_S)
{
	auto interval = FindInterval(l, lStep, _x, _S);
	auto l1 = interval.first + 0.381966011*(interval.second - interval.first);
	auto l2 = interval.second - 0.381966011*(interval.second - interval.first);
	auto f_l1 = f((_x + _S*l1).vec);
	auto f_l2 = f((_x + _S*l2).vec);

	for (auto i = 0; i < M1 && abs(interval.second - interval.first) > eps1; i++)
		if (f_l1 < f_l2)
		{
			interval.second = l2;
			l2 = l1;
			f_l2 = f_l1;
			l1 = interval.first + 0.381966011*(interval.second - interval.first);
			f_l1 = f((_x + _S*l1).vec);
		}
		else
		{
			interval.first = l1;
			l1 = l2;
			f_l1 = f_l2;
			l2 = interval.second - 0.381966011*(interval.second - interval.first);
			f_l2 = f((_x + _S*l2).vec);
		}

	return (interval.first + interval.second) / 2;
}

