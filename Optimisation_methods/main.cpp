#include <ctime>

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
#include "AdaptiveCasualMethod.h"

void test(func f, vector<func> G, Solver *s, string outFileName)
{
	double r = 1e12;
	double C = 10.;
	double penalty_eps = 1e-16;
	int M = 100;
	ofstream fout(outFileName);
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis1(-10, -1);
	uniform_real_distribution<> dis2(-10, -3);
	for (int i = 0; i <= 8; i++)
	//for (r = 1e-20; r <= 1e+20; r*=100)
	{
		Vertex x(vector<double>{dis1(gen), dis2(gen)});
		PenaltyMethod pm(C, r, penalty_eps, M, x, f, G, s);
		auto res = pm.Calc();
		fout << x << "\t" << res << "\t" << f(res.vec) << "\t" << s->GetFuncCnt() << "\t" << r  << endl;
	}
	//PenaltyMethod pm(C, r, penalty_eps, M, x, f, G, s);
	//auto res = pm.Calc();
	//fout << res << "\t" << f(res.vec) <<  "\t" << s->GetFuncCnt() << "\t" << r << endl;
}

void test2(func f, vector<func> G, Solver *s, string outFileName)
{
	Vertex x(vector<double>{0.2, 0.2, 0.2, 0.2, 0.2});
	//Vertex x(vector<double>{0.2, 0.2});
	double r = 1e-4;
	double C = 10.;
	double penalty_eps = 1e-13;
	int M = 100;
	Vertex L(vector<double>{2, 0.1, 0.1, 0.1, 0});
	Vertex R(vector<double>{100, 2, 0.45, 0.45, 1});
	//Vertex L(vector<double>{-512, -512});
	//Vertex R(vector<double>{512, 512});
	ofstream fout(outFileName);
	fout.setf(ios::scientific);
	fout.precision(16);
	//for (C = 1e-40; C <= 10; C *= 10)

	for (int i = 100000; i <= 500000; i += 100000)
	{
		for (size_t j = 0; j < 6; j++)
		{
			PenaltyMethod pm(C, r, penalty_eps, M, x, f, G, s);
			GlobalMinimizator gm(L, R);

			auto res = gm.Calc(&pm, f, i);

			cout << j << endl;
			fout << res << "\t" << f(res.vec) << "\t" << gm.GetFuncCount() << "\t" << i << endl;
		}
		cout << i << endl;
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
	//	[](vector<double> vec) { return vec[0] + 1; },
	//	[](vector<double> vec) { return vec[1] + 3; },
	//	[](vector<double> vec) { return -10 - vec[0]; },
	//	[](vector<double> vec) { return -10 - vec[1]; }
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

	//cout.setf(ios::scientific);
	//cout.precision(16);
	//cout << g(Vertex(vector<double>{
	//		2.000003699151131e+00,
	//		1.000000000000000e-01,
	//		1.000023521497010e-01,
	//		4.472463285653688e-01,
	//		1.000000000000000e+00
	//}).vec) << endl;
	//system("pause");

	ifstream fin("input/input.txt");


	//Solver *BoksSolver = new BoksMethod(fin);
	//auto x = BoksSolver->Calc(f1, Vertex(vector<double>{50, 1, 0.3, 0.3, 0.50}));
	//cout << f(x.vec) << endl << x << endl << BoksSolver->GetFuncCnt() << endl;
	//system("pause");


//	test2(f, G, new NonlinearConjugateGradientMethod(fin), "output/f_cgm_output.txt");
	//cout << "f NonlinearConjugateGradientMethod" << endl;

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f, G, new NewtonMethod(fin), "output/f_newton_output.txt");
	//cout << "f NewtonMethod" << endl;

	//fin.clear();
	//fin.seekg(0, ios::beg);
//	test(f1, G, new StaticalGradMethod(fin), "output/StaticalGrad_output.txt");
	//test
	//cout << "f StaticalGradMethod" << endl;

	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f, G, new RosenbrockMethod(fin), "output/f_rosenbrock_output.txt");
	//cout << "f RosenbrockMethod" << endl;

	////fin.clear();
	////fin.seekg(0, ios::beg);
	//test2(g, G, new NewtonMethod(fin), "output/g_NewtonMethod_output.txt");
	//fin.clear();
	//fin.seekg(0, ios::beg);
	//test2(f, G, new NewtonMethod(fin), "output/f_NewtonMethod_output.txt");
	////cout << "f AdaptiveCasualMethod" << endl;

	//fin.clear();
	////fin.seekg(0, ios::beg);
	test2(g, G, new AdaptiveCasualMethod(fin), "output/AdaptiveCasualMethod_output.txt");
	fin.clear();
	fin.seekg(0, ios::beg);
	test2(f, G, new AdaptiveCasualMethod(fin), "output/AdaptiveCasualMethod1_output.txt");
	//fin.clear();
	//fin.seekg(0, ios::beg);
//	test(f2, G, new RosenbrockMethod(fin), "output/AdaptiveCasualMethod2_output.txt");

	//test2(f, G, new LevenbergMarquardtMethod(fin), "output/f_LM_output.txt");
	//cout << "f LevenbergMarquardtMethod" << endl;
}
