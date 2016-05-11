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
	auto prevF = 0., nextF = 10.;
	for (auto i = 0; i < M && (grad.norm() > eps || abs(nextF - prevF) < 1e-6); i++)
	{
		prevX = x;
		prevF = f(prevX.vec);
		++funcCnt;
		bool flag;
		Hessian();
		auto Hf = H;
		do
		{
			WeightH(_lambda);
			Inversion();
			x = prevX + grad*H;
			nextF = f(x.vec);
			//лямбда вылетает в бесконечность в случаях с малой разницей функций
			flag = nextF < prevF || abs(nextF - prevF) < 1e-6;
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
