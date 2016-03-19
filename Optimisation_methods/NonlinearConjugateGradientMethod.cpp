#include "NonlinearConjugateGradientMethod.h"



NonlinearConjugateGradientMethod::NonlinearConjugateGradientMethod()
{
}


NonlinearConjugateGradientMethod::~NonlinearConjugateGradientMethod()
{
}

void NonlinearConjugateGradientMethod::NCGM(func _f, Vertex _x)
{
	_gss = GoldenSectionSearch(_f, eps);
	f = _f; x = _x; N = x.vec.size();
	Grad();
	S = grad;
	while (S.norm() > eps)
	{
		lambda = _gss.GSS(x, S);
		x = x + S*lambda;
		auto gn = grad.norm();
		auto w1 = gn*gn;
		Grad();
		gn = grad.norm();
		w = gn*gn / w1;
		S = grad + S*w;
	}
}

void NonlinearConjugateGradientMethod::Grad()
{
	for (auto i = 0; i < N; i++)
	{
		x1 = x; x1.vec[i] += h;
		grad.vec[i] = -(f(x1.vec) - f(x.vec)) / h;
	}
}
