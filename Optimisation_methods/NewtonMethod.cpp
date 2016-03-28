#include "NewtonMethod.h"

NewtonMethod::NewtonMethod() : eps(1e-15), h(1), h1(1e-8)
{
}

NewtonMethod::~NewtonMethod()
{
}

Vertex NewtonMethod::NM(func _f, Vertex &_x)
{
	f = _f; x = _x; N = x.vec.size();
	Init();
	Grad();
	while (grad.norm() > eps) //разобраться с минусами
	{
		Hessian();
		Inversion();
		if (_eigenvalues.FindEigenvalues(H1))
		{
			d = grad*H;
			lambda = _gss.GSS(x, d);
			x = x - d*lambda;
		}
		else
		{
			lambda = _gss.GSS(x, grad);
			x = x - grad*lambda;
		}
		Grad();
	}
	return x;
}


void NewtonMethod::Grad()
{
	auto fx = f(x.vec);
	for (auto i = 0; i < N; i++)
	{
		x1 = x; x1.vec[i] += h1;
		grad.vec[i] = (f(x1.vec) - fx) / h1;
	}
}

void NewtonMethod::Hessian()
{
	auto fx = f(x.vec);
	for (auto i = 0; i < N; i++)
	{
		for (auto j = 0; j < N; j++)
		{
			x1 = x2 = x3 = x;
			x1.vec[i] += h;		x1.vec[j] += h;
			x2.vec[i] += h;		x3.vec[j] += h;
			H[i][j] = (f(x1.vec) - f(x2.vec) - f(x3.vec) + fx) / (h*h);
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

void NewtonMethod::Init()
{
	_gss = GoldenSectionSearch(f, eps);
	x1 = Vertex(N);
	x2 = Vertex(N);
	x3 = Vertex(N);
	grad = Vertex(N);
	d = Vertex(N);
	H = vector<vector<double>>(N, vector<double>(N));
	H1 = vector<vector<double>>(N, vector<double>(N));
}
