#include "Solver.h"
#include "GoldenSectionSearch.h"

enum mod{ PolakRibiere , FletcherReeves};

int main()
{
	Solver s;
	auto _v = vector<double>(2);
	Vertex v(_v);
	auto f = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0])
		+ 100 * (vec[1] - vec[0] * vec[0])*(vec[1] - vec[0] * vec[0]); };
	auto f2 = [](vector<double> vec) { return 2 * vec[0] * vec[0] + 2 * vec[1] * vec[1] + 2 * vec[0] * vec[1] + 20 * vec[0] + 10 * vec[1] + 10; };
	cout << s.NCGM(f2, v) << endl;
	system("pause");
}