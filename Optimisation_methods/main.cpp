#include "NewtonMethod.h"
#include "NonlinearConjugateGradientMethod.h"
#include "RosenbrockMethod.h"
#include "BoksMethod.h"
#include "LevenbergMarquardtMethod.h"
#include "DavidonFletcherPowellMethod.h"
#include "ThirdPearsonMethod.h"
#include "GreenshtadtMethod.h"
#include "GoldfarbMethod.h"

int main()
{
	auto f = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0])
		+ 100 * (vec[1] - vec[0] * vec[0])*(vec[1] - vec[0] * vec[0]); };


	auto f2 = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0])
		+ 5 * (vec[1] - vec[0])*(vec[1] - vec[0]); };

	auto f3 = [](vector<double> vec) { return vec[0] * vec[0] + vec[1] * vec[1]; };

	ifstream fin("input.txt");
	Solver *s = new GoldfarbMethod(fin);

	auto res = s->Calc(f3);
	cout << res << endl << "f = " << f3(res.vec) << endl << '\a';
	system("pause");
}
