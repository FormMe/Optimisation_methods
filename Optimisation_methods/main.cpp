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
#include "GlobalMinimizator.h"
#include "StaticalGradMethod.h"
#include <ctime>

void test(func f, vector<func> G, Solver *s, string outFileName)
{
	ifstream fin2("input/input_penalty.txt");
	PenaltyMethod pm(fin2, f, G, s);
	auto res = pm.Calc();
	ofstream fout(outFileName);
	fout << "f" << res << "\t = " << f(res.vec) << "\t" << s->GetFuncCnt() << endl;
}

void test2(func f, vector<func> G, Solver *s, string outFileName)
{
	Vertex x(vector<double>{0.2, 0.2, 0.2, 0.2, 0.2});
	//Vertex x(vector<double>{0.2, 0.2});
	double r = 1e+20;
	double C = 10.;
	double penalty_eps = 1e-15;
	int M = 100;
	Vertex L(vector<double>{2, 0.1, 0.1, 0.1, 0});
	Vertex R(vector<double>{100, 2, 0.45, 0.45, 1});
	//Vertex L(vector<double>{-512, -512});
	//Vertex R(vector<double>{512, 512});
	ofstream fout(outFileName);
	//for (C = 1e-40; C <= 10; C *= 10)

	for (int i = 1000000; i < 10000000; i += 1000000)
	{
		PenaltyMethod pm(C, r, penalty_eps, M, x, f, G, s);
		GlobalMinimizator gm(L, R);

		auto start_time = clock(); // начальное время

		auto res = gm.Calc(&pm, f, i);

		auto end_time = clock(); // конечное время
		auto search_time = end_time - start_time; // искомое время

		fout << i <<"\tf" << res << "\t = " << f(res.vec) << "\t" << gm.GetFuncCount() << "\t" << search_time / 1000 << endl;

	}


	//for (r = 1e-40; r <= 1e+40; r *= 10)
	//{
	//	PenaltyMethod pm(C, r, penalty_eps, M, x, f, G, s);
	//	GlobalMinimizator gm(L, R);
	//	auto res = gm.Calc(&pm, f, 100);
	//	fout << C << "\t" << r << "\t" << res << "\t" << pm.GetFuncCount() << endl;
	//}
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

	auto f5 = [](vector<double> vec)
	{
		return -(vec[1] + 47)*
			sin(sqrt(abs(vec[0] / 2 + vec[1] + 47)))
			- vec[0] * sin(sqrt(abs(vec[0] - (vec[1] + 47))));
	};

	auto f6 = [](vector<double> vec)
	{
		return abs(
			sin(vec[0])*cos(vec[1])*
			exp(abs(1 - sqrt(vec[0] * vec[0] + vec[1] * vec[1]) / _Pi)));
	};

	auto f7 = [](vector<double> vec)
	{
		return -20 * exp(-0.2*sqrt(0.5*(vec[0] * vec[0] + vec[1] * vec[1])))
			- exp(0.5*(cos(2 * _Pi*vec[0]) + cos(2 * _Pi*vec[1]))) + _Exp1 + 20;
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
		return -(vec[0] * vec[1] * vec[4]) / (1 - vec[2]) -
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
			(1 - vec[4])*(1 + vec[3]) / (1 - vec[3]);
	};


	//auto G = vector<func>({
	//	[](vector<double> vec) { return vec[0] -512; },
	//	[](vector<double> vec) { return vec[1] -512; },
	//	[](vector<double> vec) { return -512 - vec[0]; },
	//	[](vector<double> vec) { return -512 - vec[1]; }
	//});


	auto G = vector<func>({
		[](vector<double> vec) { return vec[0] - 100.0; },
		[](vector<double> vec) { return vec[1] - 2.0; },
		[](vector<double> vec) { return vec[2] - 0.45; },
		[](vector<double> vec) { return vec[3] - 0.45; },
		[](vector<double> vec) { return vec[4] - 1.0; },
		[](vector<double> vec) { return 2.0 - vec[0]; },
		[](vector<double> vec) { return 0.1 - vec[1]; },
		[](vector<double> vec) { return 0.1 - vec[2]; },
		[](vector<double> vec) { return 0.1 - vec[3]; },
		[](vector<double> vec) { return 0.0 - vec[4]; },
	});


	ifstream fin("input/input.txt");


	//Solver *BoksSolver = new BoksMethod(fin);
	//cout << BoksSolver->Calc(f1, Vertex(vector<double>{0, 0})) << endl;
	//cout << BoksSolver->GetFuncCnt() << endl;
	//system("pause");


	//test2(g, G, new NonlinearConjugateGradientMethod(fin), "output/cgm_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(g, G, new NewtonMethod(fin), "output/newton_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	test2(f, G, new StaticalGradMethod(fin), "output/StaticalGrad_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f, G, new RosenbrockMethod(fin), "output/rosenbrock_output.txt");

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f, G, new LevenbergMarquardtMethod(fin), "output/LM_output.txt");

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

	/* дальше с файла*/

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
