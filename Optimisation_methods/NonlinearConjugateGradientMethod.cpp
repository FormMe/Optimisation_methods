#include "NonlinearConjugateGradientMethod.h"



NonlinearConjugateGradientMethod::NonlinearConjugateGradientMethod()
{
}


NonlinearConjugateGradientMethod::~NonlinearConjugateGradientMethod()
{
}

Vertex NonlinearConjugateGradientMethod::NCGM(func _f, Vertex _x)
{
	eps = 1e-8; h = 1e-8;
	_gss = GoldenSectionSearch(_f, eps);
	f = _f; x = _x; N = x.vec.size();
	x = Vertex(N);
	x1 = Vertex(N);
	grad = Vertex(N);
	S = Vertex(N);
	Grad();
	S = grad;
	while (S.norm() > eps)
	{
		lambda = _gss.GSS(x, S);
		x = x + S*lambda;
		auto gn = grad.norm();
		Grad();
		auto gn1 = grad.norm();
		w = (gn1*gn1) / (gn*gn);
		S = grad + S*w;
	}
	return x;
}

void NonlinearConjugateGradientMethod::Grad()
{
	auto fx = f(x.vec);
	for (auto i = 0; i < N; i++)
	{
		x1 = x; x1.vec[i] += h;
		grad.vec[i] = -(f(x1.vec) - fx) / h;
	}
}
