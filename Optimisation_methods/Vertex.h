#pragma once
#include <cmath>
#include <vector>
#include <cassert> 
#include <functional>
#include <algorithm>
#include <iostream>

using namespace std;

typedef  function<double(vector<double>)> func;

struct Vertex
{
	Vertex() : vec(5) {}
	Vertex(int n) : vec(n) {}
	Vertex(vector<double> &_vec) : vec(_vec) {}

	vector<double> vec;

	Vertex operator+ (Vertex _vec) const
	{
		auto v(vec);
		for (int i = 0; i < v.size(); i++)
			v[i] += _vec.vec[i];
		return Vertex(v);
	}

	Vertex operator- (Vertex _vec) const
	{
		auto v(vec);
		for (int i = 0; i < v.size(); i++)
			v[i] -= _vec.vec[i];
		return Vertex(v);
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
		auto v(vec);
		for (int i = 0; i < vec.size(); i++)
			v[i] *= c;
		return Vertex(v);
	}

	double norm()
	{
		double res = 0;
		for (int i = 0; i < vec.size(); i++)
			res += vec[i] * vec[i];
		return sqrt(res);
	}
};

