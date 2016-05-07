#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream> 
#include <algorithm>
#define e 1e-18
using namespace std;
typedef vector<vector<double>> matrix;

class LU_eigenvalues
{
public:
	LU_eigenvalues(double _eps = 1e-8) : eps(_eps) {};

	bool FindEigenvalues(matrix &M);

	int N;
	double eps;
	matrix A1;
	matrix A2;

private:
	void Decomposition();
	void Mult();
};

