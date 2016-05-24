#include "NonlinearConjugateGradientMethod.h"


NonlinearConjugateGradientMethod::NonlinearConjugateGradientMethod(ifstream& fin) :
	Solver(fin),
	grad1(Vertex(N)),
	S(Vertex(N)) {}

Vertex NonlinearConjugateGradientMethod::Calc(func _f, const Vertex &_x)
{
	x = _x;
	f = _f;
	Grad();
	S = grad;
	auto Snorm = 1.;
	for (auto i = 0; i < M &&  Snorm > eps; i++)
	{
		lambda = GSS(x, S);
		x = x + S*lambda;
		CorrectVertex(x);
		PolakRibiere();
		S = grad + S*w;
		Snorm = S.norm();
	}
	return x;
}

void NonlinearConjugateGradientMethod::PolakRibiere()
{
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
}
