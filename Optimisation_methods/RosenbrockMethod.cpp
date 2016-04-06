#include "RosenbrockMethod.h"


Vertex RosenbrockMethod::RM(func _f,  Vertex &v)
{
	N = v.vec.size();
	S = vector<Vertex>(N, Vertex(N));
	x = Vertex(N);
	//x = v;
	for (int i = 0; i < N; i++)
		S[i].vec[i] = 1;


	return {};
}
