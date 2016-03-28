#include "Solver.h"
#include "GoldenSectionSearch.h"

enum mod { PolakRibiere, FletcherReeves };

int main()
{
	Solver s;
	//auto _v = vector<double>(2);
	auto _v = vector<double>{ 0.5,1 };
	Vertex v(_v);
	auto f = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0])
		+ 100 * (vec[1] - vec[0] * vec[0])*(vec[1] - vec[0] * vec[0]); };
	auto f2 = [](vector<double> vec) { return 2 * vec[0] * vec[0] + 2 * vec[1] * vec[1] + 2 * vec[0] * vec[1] + 20 * vec[0] + 10 * vec[1] + 10; };
	auto f3 = [](vector<double> vec) { return 2 * vec[0] * vec[0] + vec[1] * vec[1] + vec[0] * vec[1]; };
	cout << s.NM(f3, v) << endl;
	system("pause");
}