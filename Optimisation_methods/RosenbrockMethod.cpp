#include "RosenbrockMethod.h"


RosenbrockMethod::RosenbrockMethod()
{
}


RosenbrockMethod::~RosenbrockMethod()
{
}

void RosenbrockMethod::RM(func _f, int n, Vertex &v)
{
	N = n;
	S = vector<Vertex>(N, Vertex(N));
	x = Vertex(N);
	//x = v;
	for (int i = 0; i < N; i++)
		S[i].vec[i] = 1;

	
}