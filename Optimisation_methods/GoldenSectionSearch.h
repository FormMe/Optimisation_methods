#pragma once
#include "Vertex.h"

class GoldenSectionSearch
{
public:
	GoldenSectionSearch();
	GoldenSectionSearch(func _f, double e) : f(_f) , eps(e) {};
	~GoldenSectionSearch();

	double GSS( Vertex &_x, Vertex &_S);

private:
	Vertex x;
	Vertex S;
	func f;
	double eps;

	pair<double, double> FindInterval(double x0, double d);
};

