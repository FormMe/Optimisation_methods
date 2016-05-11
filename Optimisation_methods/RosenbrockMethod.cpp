#include "RosenbrockMethod.h"

RosenbrockMethod::RosenbrockMethod(ifstream& fin) :
	Solver(fin),
	S(vector<Vertex>(N, Vertex(N))),
	A(vector<Vertex>(N, Vertex(N))),
	x1(Vertex(N)),
	lambda(vector<double>(N)) {}

Vertex RosenbrockMethod::Calc(func _f, const Vertex &_x)
{
	x = _x;
	f = _f;
	auto f = false;
	for (auto i = 0; i < N; i++)	
		S[i].vec[i] = 1;
	for (auto i = 0; i < M && !f; i++)
	{
		prevX = x;
		MinDirections();
		FindDirectios();
		GramSchmidtProcess();
		f = true;
		for (auto j = 0; j < N && f; j++)
			f = fabs(x.vec[0] - prevX.vec[0]) <= eps;
	}
	return x;
}

void RosenbrockMethod::MinDirections()
{
	for (auto j = 0; j < N; j++)
	{
		//слишком маламя лямбда не влияет на х
		lambda[j] = GSS(x, S[j]);
		x = x + S[j] * lambda[j];
		CorrectVertex(x);
	}
}

void RosenbrockMethod::FindDirectios()
{
	auto ind = vector<int>(N);
	for (auto k = 0; k < N; k++) ind[k] = k;
	sort(ind.begin(), ind.end(), [&](const int l, const int r) {return abs(lambda[l]) > abs(lambda[r]); });

	for (auto j = 0; j < N; ++j)
	{
		auto indx = ind[N - 1];
		A[j] = S[indx] * lambda[indx];
		for (auto m = j; m < N - 1; ++m)
		{
			indx = ind[m];
			A[j] = A[j] + S[indx] * lambda[indx];
		}
	}
}

void RosenbrockMethod::GramSchmidtProcess()
{
	S[0] = A[0] / A[0].norm();
	for (auto l = 1; l < N; l++)
	{
		auto B = A[l];
		for (auto m = 0; m < l; m++)
		{
			B = B - S[m] * (A[l] * S[m]);
		}
		auto norm = B.norm();
		if (norm > eps)
			S[l] = B / norm;
	}
}

void RosenbrockMethod::PalmerProcess()
{
	S[0] = A[0] / A[0].norm();
	for (auto l = 1; l < N; l++)
	{
		auto norm = A[l].norm();
		auto norm1 = A[l - 1].norm();
		auto z = norm*norm1*sqrt(norm1*norm1 - norm*norm);
		if (z > eps*eps)
			S[l] = (A[l] * norm1*norm1 - A[l - 1] * norm*norm) / z;
	}
}
