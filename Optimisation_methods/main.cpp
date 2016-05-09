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

void test(func f, vector<func> G, Solver *s, string outFileName)
{
	ifstream fin2("input_penalty.txt");
	PenaltyMethod pm(fin2, f, G, s);
	auto res = pm.Calc();
	ofstream fout(outFileName);
	fout << res << endl << "f = " << f(res.vec) << "\t" << s->GetFuncCnt() << '\a';
}


int main()
{
	auto f1 = [](vector<double> vec)
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

	//x - vec[0]
	//y - vec[1]
	//z - vec[2]
	//u - vec[3]
	//v - vec[4]

	auto p = [](vector<double> vec)
	{
		return (vec[4] * vec[2]) / (1 - vec[2]) +
			vec[3] * (1 - vec[4]) / (1 - vec[3]);
	};
	auto s = [](vector<double> vec)
	{
		return (vec[0] * vec[1] * vec[4]) / (1 - vec[2]) +
			(1 - vec[4]) / (1 - vec[3]);
	};
	auto t = [](vector<double> vec)
	{
		return (vec[0] * vec[4]) / (1 - vec[2] * vec[2]) +
			(1 - vec[4]) / (1 - vec[3] * vec[3]);
	};
	auto q = [](vector<double> vec)
	{
		return (vec[0] * vec[2] * vec[4]) / (1 - vec[2] * vec[2]) +
			vec[3] * (1 - vec[4]) / (1 - vec[3] * vec[3]);
	};
	auto g = [&](vector<double> vec)
	{
		return -s(vec) / (t(vec) + q(vec));
	};
	auto f = [&](vector<double> vec)
	{
		return (2 * p(vec)*s(vec)) / (t(vec) + q(vec)) +
			(vec[4] * vec[1] * (1 + vec[2])) / (1 - vec[2]) +
			(1 - vec[4])*(1 + vec[3])*(1 - vec[3]);
	};


	auto G = vector<func>({
		[](vector<double> vec) { return -1 + vec[0] + vec[1]; },
		/*	[](vector<double> vec) { return -vec[0] - 2; },
			[](vector<double> vec) { return vec[1] - 2; },
			[](vector<double> vec) { return -vec[1] - 2; }*/
	});


	//auto G = vector<func>({
	//	[](vector<double> vec) { return vec[0] - 100.0; },
	//	[](vector<double> vec) { return 2.0 - vec[0]; },
	//	[](vector<double> vec) { return vec[1] - 2.0; },
	//	[](vector<double> vec) { return 0.1 - vec[1]; },
	//	[](vector<double> vec) { return vec[2] - 0.45; },
	//	[](vector<double> vec) { return 0.1 - vec[2]; },
	//	[](vector<double> vec) { return vec[3] - 0.45; },
	//	[](vector<double> vec) { return 0.1 - vec[3]; },
	//	[](vector<double> vec) { return vec[2] - 1.0; },
	//	[](vector<double> vec) { return 0.0 - vec[2]; },
	//});


	ifstream fin("input.txt");

	//Solver *BoksSolver = new BoksMethod(fin);

	test(f1, G, new NewtonMethod(fin), "newton_output.txt");

	test(f1, G, new NonlinearConjugateGradientMethod(fin), "cgm_output.txt");

	test(f1, G, new RosenbrockMethod(fin), "rosenbrock_output.txt");

	test(f1, G, new LevenbergMarquardtMethod(fin), "LM_output.txt");

	test(f1, G, new DavidonFletcherPowellMethod(fin), "DFP_output.txt");

	test(f1, G, new ThirdPearsonMethod(fin), "third_P_output.txt");

	test(f1, G, new GreenshtadtMethod(fin), "greenshtadt_output.txt");

	test(f1, G, new GoldfarbMethod(fin), "Goldfrab_output.txt");

	system("pause");
}
