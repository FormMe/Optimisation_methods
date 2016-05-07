#include "NewtonMethod.h"
#include "NonlinearConjugateGradientMethod.h"
#include "RosenbrockMethod.h"
#include "BoksMethod.h"
#include "LevenbergMarquardtMethod.h"
#include "DavidonFletcherPowellMethod.h"
#include "ThirdPearsonMethod.h"
#include "GreenshtadtMethod.h"
#include "GoldfarbMethod.h"
#include "PenaltyMethod.h"

int main()
{
	auto f = [](vector<double> vec)
	{
		return (1 - vec[0])*(1 - vec[0])
			+ 100 * (vec[1] - vec[0] * vec[0])*(vec[1] - vec[0] * vec[0]);
	};


	auto f2 = [](vector<double> vec)
	{
		return (1 - vec[0])*(1 - vec[0])
			+ 5 * (vec[1] - vec[0])*(vec[1] - vec[0]);
	};

	auto f3 = [](vector<double> vec) { return vec[0] * vec[0] + vec[1] * vec[1]; };

	auto f4 = [](vector<double> vec)
	{
		int A1 = 1, A2 = 3;
		int a1 = 2, a2 = 1;
		int b1 = 3, b2 = 1;
		int c1 = 2, c2 = 1;
		int d1 = 3, d2 = 2;

		return -((A1 / (1 + pow(((vec[0] - a1) / b1), 2) + pow(((vec[1] - c1) / d1), 2))) 
			+ (A2 / (1 + pow(((vec[0] - a2) / b2), 2) + pow(((vec[1] - c2) / d2), 2))));

	};


	auto g = vector<func>({
		[](vector<double> vec) { return -1 + vec[0] + vec[1]; },
		/*	[](vector<double> vec) { return -vec[0] - 2; },
			[](vector<double> vec) { return vec[1] - 2; },
			[](vector<double> vec) { return -vec[1] - 2; }*/
	});

	ifstream fin("input.txt");
	ifstream fin2("input_penalty.txt");
	Solver *s = new RosenbrockMethod(fin);

	//auto x = Vertex(vector<double>{1.02535, 1.09424});
	//auto res = s->Calc(f4, x);

	PenaltyMethod pm(fin2, f4, g, s);
	auto res = pm.Calc();

	cout << res << endl << "f = " << f4(res.vec) << endl << '\a';
	system("pause");
}
