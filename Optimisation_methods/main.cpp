#include "Solver.h"

int main()
{
	Solver s;
	auto _v = vector<double>(2);
	Vertex v(_v);
	auto f = [](vector<double> vec) { return (1 - vec[0])*(1 - vec[0]) 
												+ 100 * (vec[1] - vec[0] * vec[0])*(vec[1] - vec[0] * vec[0]); };
	s.DSM(f, v, 1, 0.5, 2);

}
