#include "NewtonMethod.h"



NewtonMethod::NewtonMethod()
{
}


NewtonMethod::~NewtonMethod()
{
}

Vertex NewtonMethod::NM(func _f, Vertex _x)
{
	eps = 1e-8; h = 1;
	_gss = GoldenSectionSearch(_f, eps);
	f = _f; x = _x; N = x.vec.size();
	x1 = Vertex(N);
	x2 = Vertex(N);
	x3 = Vertex(N);
	x4 = Vertex(N);
	grad = Vertex(N);
	d = Vertex(N);
	H = vector<vector<double>>(N, vector<double>(N));
	H1 = vector<vector<double>>(N, vector<double>(N));
	Grad();
	while (grad.norm() > eps)
	{
		Hessian();
		Inversion();
		if (_eigenvalues.FindEigenvalues(H1))
		{
			d = grad*H;
			lambda = _gss.GSS(x, d);
			x = x + d*lambda;
		}
		else
		{
			lambda = _gss.GSS(x, grad);
			x = x + d*lambda;
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
		x1 = x; x1.vec[i] += h;
		grad.vec[i] = -(f(x1.vec) - fx) / h;
	}
}

void NewtonMethod::Hessian()
{
	auto dh = 4 * h*h;
	for (auto i = 0; i < N; i++)
	{
		for (auto j = 0; j < N; j++)
		{
			x1 = x2 = x3 = x4 = x;
			x1.vec[i] += h;		x1.vec[j] += h;
			x2.vec[i] += h;		x2.vec[j] -= h;
			x3.vec[i] -= h;		x3.vec[j] += h;
			x4.vec[i] -= h;		x4.vec[j] -= h;
			H[i][j] = (f(x1.vec) - f(x2.vec) - f(x3.vec) + f(x4.vec)) / dh;
			//return (f(x + hx, y + hy) - f(x + hx, y - hy) - f(x - hx, y + hy) + f(x - hx, y - hy)) / (4 * hx*hy);
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
