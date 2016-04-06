#include "Solver.h"
#include "NewtonMethod.h"

int main()
{
	Solver s = NewtonMethod();
	//auto _v = vector<double>(2);
	auto _v = vector<double>{ 0, 0 };
	Vertex v(_v);
	auto f = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0])
		+ 100 * (vec[1] - vec[0] * vec[0])*(vec[1] - vec[0] * vec[0]); };
	auto f2 = [](vector<double> vec) { return 2 * vec[0] * vec[0] + 2 * vec[1] * vec[1] + 2 * vec[0] * vec[1] + 20 * vec[0] + 10 * vec[1] + 10; };
	auto f3 = [](vector<double> vec) { return 2 * vec[0] * vec[0] + vec[1] * vec[1] + vec[0] * vec[1]; };

	cout << s.Calc(f, v) << endl;
	system("pause");
}
