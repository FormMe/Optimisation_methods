#include "NewtonMethod.h"
#include  "NonlinearConjugateGradientMethod.h"
#include "RosenbrockMethod.h"

int main()
{
	auto f = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0])
		+ 100 * (vec[1] - vec[0] * vec[0])*(vec[1] - vec[0] * vec[0]); };


	auto f2 = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0])
		+ 5 * (vec[1] - vec[0])*(vec[1] - vec[0]); };


	NewtonMethod newton("input.txt");
	NonlinearConjugateGradientMethod cgm("input.txt");
	RosenbrockMethod ros("input.txt");

	Vertex v(vector<double>{ -10, 10});
	auto res = ros.Calc(f, v);
	cout << res << endl << "f = " << f(res.vec) << endl;
	system("pause");
}
