#pragma once
#include <cmath>
#include <vector>
#include <cassert> 
#include <functional>
#include <algorithm>

using namespace std; 

typedef  function<double(vector<double>)> func;

struct Vertex
{
	Vertex()
	{
		vec = vector<double>(5);
	}
	Vertex(int n)
	{
		vec = vector<double>(n);
	}
	Vertex(vector<double> &_vec)
	{
		vec = _vec;
	}
	~Vertex()
	{
		vec.~vector();
	}

	vector<double> vec;

	Vertex operator+ (Vertex _vec)
	{
		for (int i = 0; i < vec.size(); i++)
			vec[i] += _vec.vec[i];
		return Vertex(vec);
	}

	Vertex operator- (Vertex _vec)
	{
		for (int i = 0; i < vec.size(); i++)
			vec[i] -= _vec.vec[i];
		return Vertex(vec);
	}

	double operator* (Vertex _vec)
	{
		double res = 0;
		for (int i = 0; i < vec.size(); i++)
			res += vec[i] * _vec.vec[i];
		return res;
	}

	Vertex operator* (double c)
	{
		for (int i = 0; i < vec.size(); i++)
			vec[i] *= c;
		return Vertex(vec);
	}

	double norm()
	{
		double res = 0;
		for (int i = 0; i < vec.size(); i++)
			res += vec[i] * vec[i];
		return sqrt(res);
	}
};

