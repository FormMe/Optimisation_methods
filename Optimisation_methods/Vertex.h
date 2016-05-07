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
	Vertex() {}
	Vertex(int n) : vec(n) {}
	Vertex(const vector<double> &_vec) : vec(_vec) {}

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

	vector<vector<double>> TransMult(Vertex &_vec)
	{
		auto n = vec.size();
		vector<vector<double>> res(n, vector<double>(n));
		for (auto i = 0; i < n; i++)
			for (auto j = 0; j < n; j++)
				res[i][j] = vec[i] * _vec.vec[j];
		return res;
	}

	Vertex operator* (double c)
	{
		auto _v(vec);
		for (auto &v : _v)
			v *= c;
		return Vertex(_v);
	}

	Vertex operator/ (double c)
	{
		auto _v(vec);
		for (auto &v : _v)
			v /= c;
		return Vertex(_v);
	}

	Vertex operator* (vector<vector<double>> &M)
	{
		auto N = vec.size();
		auto res = Vertex(N);
		for (auto i = 0; i < N; ++i)
			for (auto j = 0; j < N; j++)
				res.vec[i] += M[i][j] * vec[j];
		return res;
	}

	Vertex TransMult(vector<vector<double>> &M)
	{
		auto N = vec.size();
		auto res = Vertex(N);
		for (auto i = 0; i < N; ++i)
			for (auto j = 0; j < N; j++)
				res.vec[i] += vec[j] * M[j][i];
		return res;
	}

	double norm()
	{
		double res = 0;
		for (auto v : vec)
			res += v*v;
		return sqrt(res);
	}

	double len(Vertex &ver)
	{
		double l = 0;
		for (auto i = 0; i < vec.size(); i++)
			l += (ver.vec[i] - vec[i])*(ver.vec[i] - vec[i]);
		return sqrt(l);
	}

	friend ostream& operator << (ostream& ostream_, const Vertex& v)
	{
		ostream_.setf(ios::scientific);
		ostream_.precision(16);
		ostream_ << "[ ";
		for (auto x : v.vec)
			ostream_ << x << ", ";
		ostream_ << " ]";
		return ostream_;
	}
};

