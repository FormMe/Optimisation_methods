#include "NewtonMethod.h"

NewtonMethod::NewtonMethod(ifstream& fin) :
	Solver(fin),
	x2(Vertex(N)),
	x3(Vertex(N)),
	H(vector<vector<double>>(N, vector<double>(N))),
	H1(vector<vector<double>>(N, vector<double>(N))) {}

Vertex NewtonMethod::Calc(func _f, const Vertex &_x)
{
	x = _x;
	f = _f;
	CorrectVertex(x);
	Grad();
	for (auto i = 0; i < M && grad.norm() > eps; i++)
	{
		Hessian();
		Inversion();
		if (_eigenvalues.FindEigenvalues(H1))
		{
			S = grad*H;
			lambda = Fibbonachi(x, S);
			x = x + S*lambda;			
		}
		else
		{
			lambda = Fibbonachi(x, grad);
			x = x + grad*lambda;
		}
		CorrectVertex(x);
		Grad();
	}
	return x;
}



void NewtonMethod::Hessian()
{
	auto fx = f(x.vec);
	++funcCnt;
	for (auto i = 0; i < N; i++)
	{
		for (auto j = i; j < N; j++)
		{
			x1 = x2 = x3 = x;
			x1.vec[i] += h1;		x1.vec[j] += h1;
			x2.vec[i] += h1;		x3.vec[j] += h1;
			auto f1 = f(x1.vec);
			auto f2 = f(x2.vec);
			auto f3 = f(x3.vec);
			H[i][j] = H[j][i] = (f1 - f2 - f3 + fx) / (h1*h1);
			funcCnt += 3;
		}
	}
}

void NewtonMethod::Inversion()
{
	double temp;
	H1.assign(H1.size(), vector<double>(N, 0));
	for (auto i = 0; i < N; i++)
		H1[i][i] = 1;
	for (auto k = 0; k < N; k++)
	{
		temp = H[k][k];
		for (auto j = 0; j < N; j++)
		{
			H[k][j] /= temp;
			H1[k][j] /= temp;
		}
		for (auto i = k + 1; i < N; i++)
		{
			temp = H[i][k];
			for (int j = 0; j < N; j++)
			{
				H[i][j] -= H[k][j] * temp;
				H1[i][j] -= H1[k][j] * temp;
			}
		}
	}
	for (auto k = N - 1; k > 0; k--)
		for (auto i = k - 1; i >= 0; i--)
		{
			temp = H[i][k];
			for (auto j = 0; j < N; j++)
			{
				H[i][j] -= H[k][j] * temp;
				H1[i][j] -= H1[k][j] * temp;
			}
		}
	H = H1;
}

