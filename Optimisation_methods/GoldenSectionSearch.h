#pragma once
#include "Vertex.h"

class GoldenSectionSearch
{
public:
	GoldenSectionSearch();
	~GoldenSectionSearch();

	double GSS(func _f, Vertex &_x, Vertex &_S);

private:
	Vertex x;
	Vertex S;
	func f;

	pair<double, double> FindInterval(double x0, double d) ;
};

