#include "LevenbergMarquardtMethod.h"

LevenbergMarquardtMethod::LevenbergMarquardtMethod(ifstream& fin) : NewtonMethod(fin) { fin >> lambda; }

LevenbergMarquardtMethod::~LevenbergMarquardtMethod()
{
}

Vertex LevenbergMarquardtMethod::Calc(func _f)
{
	f = _f;
	Grad();
	for (auto i = 0; i < M && grad.norm() > eps; i++)
	{
		auto prevX = x;
		auto prevF = f(prevX.vec);
		bool flag;
		Hessian();
		auto Hf = H;
		do
		{
			WeightH();
			Inversion();
			x = prevX + grad*H;
			flag = f(x.vec) < prevF;
			if (flag) lambda /= 2;
			else
			{
				lambda *= 2;
				H = Hf;
			}

		} while (!flag);
		Grad();
	}
	return x;
}

void LevenbergMarquardtMethod::WeightH()
{
	for (auto i = 0; i < N; ++i)
		H[i][i] += lambda;
}
