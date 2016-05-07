#include "LevenbergMarquardtMethod.h"

LevenbergMarquardtMethod::LevenbergMarquardtMethod(ifstream& fin) : NewtonMethod(fin) { fin >> lambda; }

LevenbergMarquardtMethod::~LevenbergMarquardtMethod()
{
}

Vertex LevenbergMarquardtMethod::Calc(func _f, const Vertex &_x)
{
	x = _x;
	f = _f;
	Grad();
	for (auto i = 0; i < M && grad.norm() > eps; i++)
	{
		prevX = x;
		auto prevF = f(prevX.vec);
		++funcCnt;
		bool flag;
		Hessian();
		auto Hf = H;
		do
		{
			WeightH();
			Inversion();
			x = prevX + grad*H;
			flag = f(x.vec) < prevF;
			++funcCnt;
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
