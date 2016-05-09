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
	auto _lambda = lambda;
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
			WeightH(_lambda);
			Inversion();
			x = prevX + grad*H;
			flag = abs(f(x.vec) - prevF) < eps;
			++funcCnt;
			if (flag) _lambda /= 2;
			else
			{
				_lambda *= 2;
				H = Hf;
			}

		} while (!flag);
		Grad();
	}
	return x;
}

void LevenbergMarquardtMethod::WeightH(double _lambda)
{
	for (auto i = 0; i < N; ++i)
		H[i][i] += _lambda;
}
