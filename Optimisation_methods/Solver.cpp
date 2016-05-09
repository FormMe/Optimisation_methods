#include "Solver.h"

Solver::Solver(ifstream &fin) : funcCnt(0)
{
	fin >> eps >> eps1 >> l >> lStep >> h >> h1 >> M >> M1 >> N;
	x = Vertex(N);
	prevX = Vertex(N);
	S = Vertex(N);
}


pair<double, double> Solver::FindInterval(double lambda0, double d, Vertex &_x, Vertex &_S)
{
	auto lambda1 = lambda0, lambda2 = lambda1 + d;
	funcCnt += 2;

	if (f((_x + _S*lambda0).vec) < f((_x + _S*(lambda0 + d)).vec))
		d = -d;


	funcCnt += 2;
	while (f((_x + _S*lambda1).vec) > f((_x + _S*lambda2).vec))
	{
		d *= 2;
		lambda0 = lambda1;
		lambda1 = lambda2;
		lambda2 = lambda1 + d;
		funcCnt += 2;
	}

	//if (f((_x + _S*lambda0).vec) < f((_x + _S*(lambda0 + d)).vec))
	//	d = -d;

	//do
	//{
	//	lambda1 = lambda0 + d;
	//	d *= 2;
	//	lambda2 = lambda1 + d;
	//	lambda0 = lambda1;
	//	funcCnt += 2;

	//} while (f((_x + _S*lambda1).vec) > f((_x + _S*lambda2).vec));

	//lambda0 -= d / 2;
	//пофиксить вылетание за границы
	return make_pair(min(lambda0, lambda2), max(lambda0, lambda2));
}

//GSS
double Solver::GSS(Vertex &_x, Vertex &_S)
{
	auto interval = FindInterval(l, lStep, _x, _S);
	auto l1 = interval.first + 0.381966011*(interval.second - interval.first);
	auto l2 = interval.second - 0.381966011*(interval.second - interval.first);
	auto f_l1 = f((_x + _S*l1).vec);
	auto f_l2 = f((_x + _S*l2).vec);
	funcCnt += 2;

	for (auto i = 0; i < M1 && abs(interval.second - interval.first) > eps1; i++)
		if (f_l1 < f_l2)
		{
			interval.second = l2;
			l2 = l1;
			f_l2 = f_l1;
			l1 = interval.first + 0.381966011*(interval.second - interval.first);
			f_l1 = f((_x + _S*l1).vec);
			++funcCnt;
		}
		else
		{
			interval.first = l1;
			l1 = l2;
			f_l1 = f_l2;
			l2 = interval.second - 0.381966011*(interval.second - interval.first);
			f_l2 = f((_x + _S*l2).vec);
			++funcCnt;
		}

	return (interval.first + interval.second) / 2;
}

//Fibbonachi
double Solver::Fibbonachi(Vertex& _x, Vertex& _S)
{
	auto s = FindInterval(l, lStep, _x, _S);
	double f1, f2;
	auto fib_max = (s.second - s.first) /eps1;
	long long int add_fib = 0;
	vector<long long int> fibs{ 1,1 };
	auto n = 2;
	while (fib_max > add_fib)
	{
		add_fib = fibs[n - 1] + fibs[n - 2];
		fibs.push_back(add_fib);
		n++;
	}
	n = fibs.size() - 3;

	auto l1 = s.first +
		(static_cast<double>(fibs[n]) / static_cast<double>(fibs[n + 2])) * (s.second - s.first);
	auto l2 = s.first + s.second - l1;
	f1 = f((_x + _S*l1).vec);
	f2 = f((_x + _S*l2).vec);

	for (auto k = 1; k < n; k++)
	{
		if (f1 < f2)
		{
			s.second = l2;
			l2 = l1;
			f2 = f1;
			l1 = s.first + (static_cast<double>(fibs[n - k + 1]) /
				static_cast<double>(fibs[n - k + 3])) * (s.second - s.first);
			f1 = f((_x + _S*l1).vec);
		}
		else
		{
			s.first = l1;
			l1 = l2;
			f1 = f2;
			l2 = s.first + (static_cast<double>(fibs[n - k + 2]) /
				static_cast<double>(fibs[n - k + 3])) * (s.second - s.first);
			f2 = f((_x + _S*l2).vec);
		}
	}
	return (s.first + s.second) / 2;
}



int Solver::GetFuncCnt()
{
	return funcCnt;
}