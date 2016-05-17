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
	ifstream fin2("input/input_penalty.txt");
	PenaltyMethod pm(fin2, f, G, s);
	auto res = pm.Calc();
	ofstream fout(outFileName);
	fout << res << endl << "f = " << f(res.vec) << endl << s->GetFuncCnt();
}

void test2(func f, vector<func> G, Solver *s, string outFileName)
{
	Vertex x(vector<double>{0.2, 0.2});
	//double r = 1.;
	double C = 10.;
	double penalty_eps = 1e-8;
	int M = 100;

	ofstream fout(outFileName);
	//for (auto C = 1e-40; C <= 10; C *= 10)
	//ÀÒÅÍØÍ. ÌÅÒÎÄ ÀÁÍÎÂËßÈÖÀ ÕÓÉÎÂÎ
	for (auto r = 1e-40; r <= 1e+40; r *= 10)
	{
		PenaltyMethod pm(C, r, penalty_eps, M, x, f, G, s);
		auto res = pm.Calc();
		fout << C << "\t" << r << "\t" << res << "\t" << pm.GetFuncCount() << endl;
	}
}
//
//
//void test3(func f, vector<func> G, string outFileName)
//{
//	Vertex x(vector<double>{0.2, 0.2});
//	//double r = 1.;
//	double C = 10.;
//	double penalty_eps = 1e-8;
//	int M = 100;
//
//	ofstream fout(outFileName);
//	//for (auto C = 10e-20; C <= 1e20; C *= 10)
//	for (auto r = 10e-20; r <= 1e20; r *= 10)
//	{
//		PenaltyMethod pm(C, r, penalty_eps, M, x, f, G, s);
//		auto res = pm.Calc();
//		fout << C << "\t" << r << "\t" << res << "\t" << pm.GetFuncCount() << endl;
//	}
//}


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
		[](vector<double> vec) { return vec[0] - 0.5; },
		[](vector<double> vec) { return vec[1] - 0.5; },
		[](vector<double> vec) { return -1 - vec[0]; },
		[](vector<double> vec) { return -1 - vec[1]; }
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


	ifstream fin("input/input.txt");


	//Solver *BoksSolver = new BoksMethod(fin);
	//cout << BoksSolver->Calc(f1, Vertex(vector<double>{0, 0})) << endl;
	//cout << BoksSolver->GetFuncCnt() << endl;
	//system("pause");


	//test2(f1, G, new NonlinearConjugateGradientMethod(fin), "output/cgm_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f1, G, new NewtonMethod(fin), "output/newton_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	test2(f1, G, new RosenbrockMethod(fin), "output/rosenbrock_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f1, G, new LevenbergMarquardtMethod(fin), "output/LM_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f1, G, new DavidonFletcherPowellMethod(fin), "output/DFP_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f1, G, new ThirdPearsonMethod(fin), "output/third_P_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f1, G, new GreenshtadtMethod(fin), "output/greenshtadt_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f1, G, new GoldfarbMethod(fin), "output/Goldfrab_output.txt");

	/* äàëüøå ñ ôàéëà*/

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test(f1, G, new NonlinearConjugateGradientMethod(fin), "output/cgm_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test(f1, G, new NewtonMethod(fin), "output/newton_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test(f1, G, new RosenbrockMethod(fin), "output/rosenbrock_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test(f1, G, new LevenbergMarquardtMethod(fin), "output/LM_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test(f1, G, new DavidonFletcherPowellMethod(fin), "output/DFP_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test(f1, G, new ThirdPearsonMethod(fin), "output/third_P_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test(f1, G, new GreenshtadtMethod(fin), "output/greenshtadt_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test(f1, G, new GoldfarbMethod(fin), "output/Goldfrab_output.txt");
}
