#include "GlobalMinimizator.h"
#include <random>

Vertex GlobalMinimizator::Calc(PenaltyMethod* solver, func f, int M)
{
	vector<Vertex> x_vector(M, Vertex(N));
	vector<double> Fx_vector(M);

	if(true)
	{
		
	std::generate(x_vector.begin(), x_vector.end(),
		[&]()
	{
		GetRandomX();
		return x;
	});

	auto i = 0;
	std::generate(Fx_vector.begin(), Fx_vector.end(),
		[&]()
	{
		return f(x_vector[i++].vec);
	});


	auto ind = vector<int>(M);
	for (auto k = 0; k < M; k++) ind[k] = k;
	sort(ind.begin(), ind.end(), [&Fx_vector](const int l, const int r) {return Fx_vector[l] < Fx_vector[r]; });
	auto minF_ind = ind.front();

//	auto minF_ind = min_element(Fx_vector.begin(), Fx_vector.end()) - Fx_vector.begin();


	x = x_vector[minF_ind];
	funcCnt += M;

	}
	else GetRandomX();

	auto minF = f(x.vec);
	bool isNewXfound = false;
	auto stop = 0;
	//int i = 0;
	do
	{
		x = solver->Calc(x);
		++funcCnt;
		auto fx = f(x.vec);
		if (fx < minF)
		{
			minX = x;
			minF = fx;
		}

		funcCnt += solver->GetFuncCount();
	//	GetRandomX();

		isNewXfound = false;
		for (int j = 0; j < M && !isNewXfound; j++)
		{
			GetRandomX();
			isNewXfound = f(x.vec) < minF;
			++funcCnt;
		}
	} while (isNewXfound);

	return minX;
}

void GlobalMinimizator::GetRandomX()
{
	random_device rd;
	mt19937 gen(rd());
	for (auto i = 0; i < N; i++)
	{
		uniform_real_distribution<> dis(L.vec[i], R.vec[i]);
		x.vec[i] = dis(gen);
	}
}

int GlobalMinimizator::GetFuncCount()
{
	return funcCnt;
}
