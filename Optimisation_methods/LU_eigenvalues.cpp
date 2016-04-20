#include "LU_eigenvalues.h"
LU_eigenvalues::~LU_eigenvalues()
{
}


void LU_eigenvalues::Decomposition()
{
	for (int i = 0; i < N; i++)
	{
		double sumD = 0;
		for (int j = 0; j < i; j++)
		{
			double sumL = 0, sumU = 0;
			for (int m = 0; m < j; m++)
			{
				sumL += A1[i][m] * A1[m][j];
				sumU += A1[j][m] * A1[m][i];
			}
			A1[i][j] -= sumL;
			A1[i][j] /= A1[j][j];
			A1[j][i] -= sumU;
			sumD += A1[i][j] * A1[j][i];
		}
		A1[i][i] -= sumD;
	}
}

void LU_eigenvalues::Mult()
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			A2[i][j] = 0;
			for (int l = j <= i ? i : j; l < N; l++)
				if (l == j)
					A2[i][j] += A1[i][l];
				else
					A2[i][j] += A1[i][l] * A1[l][j];
		}
	A1 = A2;
}

bool LU_eigenvalues::FindEigenvalues(matrix &M)
{
	A1 = M;
	N = A1.size();
	A2.resize(N);
	for (int i = 0; i < N; i++)
		A2[i].resize(N);
	bool R = false;
	for (int k = 0;!R; k++)
	{
		try { Decomposition(); }
		catch (char* err) { cout << k << " " << err; return false; }
		Mult();
		R = true;
		for (int i = 0; i < N && R; i++)
		{
			double norm = 0;
			for (int j = i + 1; j < N; j++)
				norm += A1[j][i] * A1[j][i];
			R = !(sqrt(norm) > eps);
		}
		
	}
	R = true;
	for (int i = 0; i < N && R; i++)
		R = A1[i][i] >= 0;
	return R;
}