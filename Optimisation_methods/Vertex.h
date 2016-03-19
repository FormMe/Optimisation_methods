#pragma once
#include <cmath>
#include <vector>
#include <cassert> 
#include <functional>
#include <algorithm>
#include <iostream>

using namespace std;

typedef  function<double(vector<double>)> func;
typedef  function<double(double)> func2;

struct Vertex
{
	Vertex() : vec(5) {}
	Vertex(int n) : vec(n) {}
	Vertex(vector<double> &_vec) : vec(_vec) {}

	vector<double> vec;

	Vertex operator+ (Vertex &_vec) 
	{
		auto v(vec);
		for (auto i = 0; i < v.size(); i++)
			v[i] += _vec.vec[i];
		return Vertex(v);
	}

	Vertex operator- (Vertex &_vec) 
	{
		auto v(vec);
		for (auto i = 0; i < v.size(); i++)
			v[i] -= _vec.vec[i];
		return Vertex(v);
	}

	double operator* (Vertex &_vec)
	{
		double res = 0;
		for (auto i = 0; i < vec.size(); i++)
			res += vec[i] * _vec.vec[i];
		return res;
	}

	Vertex operator* (double c) 
	{
		auto v(vec);
		for (auto &v : vec)
			v *= c;
		return Vertex(v);
	}

	double norm()
	{
		double res = 0;
		for (auto v : vec)
			res += v*v;
		return sqrt(res);
	}
};

