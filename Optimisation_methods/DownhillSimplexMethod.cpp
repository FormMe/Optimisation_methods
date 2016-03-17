#include "DownhillSimplexMethod.h"

DownhillSimplexMethod::DownhillSimplexMethod()
{
}

DownhillSimplexMethod::~DownhillSimplexMethod()
{
}

void DownhillSimplexMethod::DSM(func _f, Vertex &v, double _alpha, double _beta, double _gamma)
{
	f = _f;
	alpha = _alpha;
	beta = _beta;
	gamma = _gamma;
	Initialization(v);
	while (!QuitCase())
	{
		MinMax();
		Average();
		Reflection();
		if (f_reflected < F[min_ind])
			Contraction();

		if ((F[min_ind] < f_reflected) && (f_reflected <= F[g_ind]))
		{	
			smplx[max_ind] = reflected;
			F[max_ind] = f_reflected;
		}

		auto flag = false;
		if ((F[max_ind] < f_reflected) && (f_reflected <= F[g_ind]))
		{	
			flag = true;
			swap(F[max_ind], f_reflected);
			swap(smplx[max_ind], reflected);
		}

		if (f_reflected > F[max_ind])	flag = true;
		if (flag) 
			Expansion();
	}
}

void DownhillSimplexMethod::Initialization(Vertex &v)
{
	N = v.vec.size();
	double t = 1,
		q = t / N / sqrt(2),
		w = sqrt(N + 1) - 1,
		d1 = q*(w + N),
		d2 = q * w;
	smplx = vector<Vertex>(N + 1, v);
	F = vector<double>(N + 1);
	max_ind = min_ind = g_ind = 0;
	contracted = Vertex(N);
	reflected = Vertex(N);
	expansed = Vertex(N);
	average = Vertex(N);
	for (auto i = 0; i < N + 1; i++)
	{
		for (auto j = 0; j < N && i != 0; j++)
			smplx[i].vec[j] = (i - 1 != j) ? d2 : d1;
		F[i] = f(smplx[i].vec);
	}
}

void DownhillSimplexMethod::MinMax()
{
	auto result = minmax_element(F.begin(), F.end());
	min_ind = result.first - F.begin();
	max_ind = result.second - F.begin();
	for (auto i = 0; i < N + 1; i++)
		if (F[g_ind] > F[min_ind] && F[g_ind] < F[max_ind])
			g_ind = i;
}

void DownhillSimplexMethod::Average()
{
	for (auto i = 0; i < N + 1; i++)
		if (i != max_ind)
			average = average + smplx[i];
	average = average*(1.0 / N);
	f_average = f(average.vec);
}


void DownhillSimplexMethod::Reflection()
{
	reflected = average*(alpha + 1) - smplx[max_ind] * alpha;
	f_reflected = f(reflected.vec);
}

void DownhillSimplexMethod::Contraction()
{
	contracted = average*(1 - gamma) + reflected*gamma;
	f_contracted = f(contracted.vec);
	if (f_contracted < F[min_ind])
	{
		smplx[max_ind] = contracted;
		F[max_ind] = f_contracted;
	}
	else
	{
		smplx[max_ind] = reflected;
		F[max_ind] = f_reflected;
	}
}

void DownhillSimplexMethod::Expansion()
{
	expansed = average*(1 - beta) + smplx[max_ind] * beta;
	f_expansed = f(expansed.vec);
	if (f_expansed < F[max_ind])
	{
		swap(F[max_ind], f_expansed);
		swap(smplx[max_ind], expansed);
	}
	else
		Reduction();
}

void DownhillSimplexMethod::Reduction()
{
	auto min = smplx[min_ind];
	for (auto i = 0; i < N + 1; i++)
	{
		smplx[i] = min + (smplx[i] - min)*0.5;
		F[i] = f(smplx[i].vec);
	}
}

bool DownhillSimplexMethod::QuitCase()
{
	double r = 0;
	for (int i = 0; i < N; i++)
		r += (F[i] - f_average)*(F[i] - f_average);
	auto res = sqrt(r / N);
	return res < eps;
}