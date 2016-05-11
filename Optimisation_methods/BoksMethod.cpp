#include "BoksMethod.h"
#include <numeric>

BoksMethod::BoksMethod(ifstream& fin) :
	Solver(fin),
	cmplx(vector<Vertex>(K)),
	F(vector<double>(K)),
	average(Vertex(N)),
	K(2*N)
{
	fin >> alpha;
	for (auto i = 0; i < N; ++i)
		fin >> L.vec[i];
	for (auto i = 0; i < N; ++i)
		fin >> R.vec[i];
}

Vertex BoksMethod::Calc(func _f)
{
	f = _f;
	InitCmplx();
	for (auto i = 0; i < M && !QuitCase(); i++)
	{
		max_ind = max_element(F.begin(), F.end()) - F.begin();
		Average();
		x = average + (average - cmplx[max_ind])*alpha;
		while (true)
		{
			CorrectVertex(x);
			auto fr = f(x.vec);
			++funcCnt;
			if (fr < F[max_ind])
			{
				cmplx[max_ind] = x;
				F[max_ind] = fr;
				break;
			}
			else x = (x + average) / 2;
		}
	}

	return cmplx[min_element(F.begin(), F.end()) - F.begin()];
}


void BoksMethod::InitCmplx()
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(0, 1);

	cmplx[0] = x;
	F[0] = f(cmplx[0].vec);
	generate(cmplx.begin() + 1, cmplx.end(), [&]() {return L + (R - L)* dis(gen); });
	auto i = 1;
	generate(F.begin() + 1, F.end(), [&]() {return f(cmplx[i++].vec); });
	funcCnt += N - 1;
}

void BoksMethod::Average()
{
	average = Vertex(N);
	for (auto i = 0; i < K; i++)
		if (i != max_ind)
			average = average + cmplx[i];
	average = average / (K - 1);

}

bool BoksMethod::QuitCase()
{
	auto ff = accumulate(F.begin(), F.end(), 0) / K;
	auto d1 = sqrt(accumulate(F.begin(), F.end(), 0,
				[&ff](double a, double b)
				{
					return (a - ff)*(a - ff) + (b - ff)*(b - ff);
				})) 
				/ (K - 1);

	auto d2 = 0.0;
	for (auto i = 0; i < K; i++)
		for (auto j = i + 1; j < K; j++)
		{
			auto len = cmplx[i].len(cmplx[j]);
			if (len > d2) d2 = len;
		}
	return d1 < eps && d2 < eps1;
}


