#include "Solver.h"

int main()
{
	Solver s;
	Vertex v(5);
	auto f = [](vector<double>) { return 1; };
	s.DSM(f, v, 1, 1, 1);
}