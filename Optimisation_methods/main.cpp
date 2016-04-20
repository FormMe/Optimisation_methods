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


	ifstream fin("input.txt");
	Solver *s = new RosenbrockMethod(fin);
	
	Vertex l(vector<double>{ 0.5, 0});
	Vertex r(vector<double>{ 1.5, 1.5});

	auto res = s->Calc(f);
	cout << res << endl << "f = " << f(res.vec) << endl << '\a';
	system("pause");
}
