#include "GlobalMinimizator.h"
#include <random>

Vertex GlobalMinimizator::Calc(PenaltyMethod* solver, func f, int M)
{
	vector<Vertex> x_vector(M, Vertex(N));
	vector<double> Fx_vector(M);
	generate(x_vector.begin(), x_vector.end(),
	[&]()
	{
		GetRandomX();
		return x;
	});
	auto i = 1;
	generate(Fx_vector.begin() + 1, Fx_vector.end(), 
	[&]()
	{
		return f(x_vector[i++].vec);
	});

	auto minF_ind = min_element(Fx_vector.begin(), Fx_vector.end()) - Fx_vector.begin();
	x = x_vector[minF_ind];

	auto isNewXfound = false;
	do
	{
		minX = solver->Calc(x);
		auto minF = f(minX.vec);
		for (size_t i = 0; i < M && !isNewXfound; i++)
		{
			GetRandomX();
			isNewXfound = f(x.vec) < minF;
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
