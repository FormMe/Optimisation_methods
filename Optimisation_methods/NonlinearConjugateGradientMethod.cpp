#include "NonlinearConjugateGradientMethod.h"


Vertex NonlinearConjugateGradientMethod::Calc(func _f, Vertex _x)
{
	f = _f; x = _x; N = x.vec.size();
	x1 = Vertex(N);
	grad = Vertex(N);
	grad1 = Vertex(N);
	Grad();
	S = grad;
	for (auto i = 0; i < M && S.norm() > eps; i++)
	{
		lambda = GSS(x, S);
		x = x + S*lambda;
		FletcherReeves();
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

void NonlinearConjugateGradientMethod::PolakRibiere()
{
	//grad1 = grad;
	//Grad();
	//w = grad*(grad - grad1) / (grad1*S);
	grad1 = grad;
	Grad();
	w = grad*(grad - grad1) / (grad1*grad1);
}

void NonlinearConjugateGradientMethod::FletcherReeves()
{
	auto gn = grad.norm();
	Grad();
	auto gn1 = grad.norm();
	w = (gn1*gn1) / (gn*gn);

	//grad1 = grad;
	//Grad();
	//w = (grad*grad) / (grad1*grad1);
}
