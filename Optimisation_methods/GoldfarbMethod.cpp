#include "GoldfarbMethod.h"



GoldfarbMethod::GoldfarbMethod(ifstream &fin) : QuasiNewtonMethod(fin)
{
}


GoldfarbMethod::~GoldfarbMethod()
{
}

void GoldfarbMethod::ApproachMatrix()
{
	auto coef = dg.TransMult(H)*dg;
	Mult(dx.TransMult(dg), H);
	auto A = (dg*H).TransMult(dx);
	auto B = Prod;
	auto coef2 = 1 + dg*dx / coef;
	Mult((dg*H).TransMult(dg), H);

	for (auto i = 0; i < N; i++)
		for (auto j = 0; j < N; j++)
			H[i][j] += (A[i][j] + B[i][j] - Prod[i][j] * coef2 / coef) / coef;
}
