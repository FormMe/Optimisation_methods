#pragma once
#include "Vertex.h"


class RosenbrockMethod
{
public:
	RosenbrockMethod();
	~RosenbrockMethod();

	func f;
	int N;
	Vertex x;
	double eps;
	vector<Vertex> S;
	double lambda;

	void RM(func _f, int n, Vertex &v);
};

