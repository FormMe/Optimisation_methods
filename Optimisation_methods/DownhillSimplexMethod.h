#pragma once
#include "Vertex.h"

class DownhillSimplexMethod
{
public:
	DownhillSimplexMethod();
	~DownhillSimplexMethod();

	void DSM(func _f, Vertex &v, double alpha, double beta, double gamma);

private:
	func f;
	int N;
	double alpha;
	double beta;
	double gamma;
	double eps;

	vector<Vertex> smplx;
	vector<double> F;
	
	int min_ind;
	Vertex min;
	double f_min;

	int max_ind;
	Vertex max;
	double f_max;

	Vertex average;
	double f_average;

	Vertex reflected;
	double f_reflected;

	Vertex contracted;
	double f_contracted;

	Vertex expansed;
	double f_expansed;

	void Initialization(Vertex &v);
	void MinMax();
	void Average();
	bool QuitCase();
	void Reflection();
	void Expansion();
	void Contraction();
	void Reduction();
};

