#include "NewtonMethod.h"
#include  "NonlinearConjugateGradientMethod.h"
#include "RosenbrockMethod.h"
#include "BoksMethod.h"

int main()
{
	auto f = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0])
		+ 100 * (vec[1] - vec[0] * vec[0])*(vec[1] - vec[0] * vec[0]); };


	auto f2 = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0])
		+ 5 * (vec[1] - vec[0])*(vec[1] - vec[0]); };


	NewtonMethod newton("input.txt");
	NonlinearConjugateGradientMethod cgm("input.txt");
	RosenbrockMethod ros("input.txt");
	BoksMethod boks("input.txt");

	Vertex v(vector<double>{ 0, 0});
	Vertex l(vector<double>{ 0.5, 0});
	Vertex r(vector<double>{ 1.5, 1.5});

	auto res = boks.Calc(f, v, l, r, 1.3);
	cout << res << endl << "f = " << f(res.vec) << endl << '\a';
	system("pause");
}
