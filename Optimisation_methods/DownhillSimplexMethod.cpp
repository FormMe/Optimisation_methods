#include "DownhillSimplexMethod.h"


DownhillSimplexMethod::DownhillSimplexMethod()
{
}


DownhillSimplexMethod::~DownhillSimplexMethod()
{
}

void DownhillSimplexMethod::DSM(func _f, Vertex &v, double _alpha, double _beta, double _gamma)
{
	N = v.vec.size();
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
		if (f_reflected < f_max)
		{
			if (f_reflected < f_min)
			{
				Contraction();
			}
			else
			{
				Expansion();
				if (f_expansed < f_max)
				{
					smplx[max_ind] = expansed;
					F[max_ind] = f_expansed;
				}
				else
					Reduction();
			}
		}
		else
			Reduction();
	}
}

void DownhillSimplexMethod::Initialization(Vertex &v)
{
	double t = 0.1,
		q = t / N / sqrt(2),
		w = sqrt(N + 1) - 1,
		d1 = q*(w + N),
		d2 = q / w;
	smplx = vector<Vertex>(N + 1,v);
	F = vector<double>(N + 1);
	max_ind = 0; min_ind = 0;
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N && i != 0; j++)
			smplx[i].vec[j] = (i - 1 != j) ? d2 : d1;
		F[i] = f(smplx[i].vec);
	}
}

void DownhillSimplexMethod::MinMax()
{
	auto result = minmax_element(F.begin(), F.end());
	min_ind = result.first - F.begin();
	max_ind = result.second - F.begin();
	min = smplx[min_ind];
	max = smplx[max_ind];
}

void DownhillSimplexMethod::Average()
{
	for (int i = 0; i < N; i++)
		if (i != max_ind)
			average = average + smplx[i];
	average = average*(1 / N);
	f_average = f(average.vec);
}

bool DownhillSimplexMethod::QuitCase()
{
	double r = 0;
	for (int i = 0; i < N; i++)
		r += (F[i] - f_average)*(F[i] - f_average);
	return sqrt(r / (N + 1)) < eps;
}

void DownhillSimplexMethod::Reflection()
{
	reflected = average*(alpha + 1) - max * alpha;
	f_reflected = f(reflected.vec);
}

void DownhillSimplexMethod::Contraction()
{
	contracted = average*(1 - gamma) + reflected*gamma;
	f_contracted = f(contracted.vec);
   if (f_contracted < f_min)
   {
      smplx[min_ind] = contracted;
      F[min_ind] = f_contracted;
   }
   else
   {
      smplx[min_ind] = reflected;
      F[min_ind] = f_reflected;
   }
}

void DownhillSimplexMethod::Expansion()
{
	expansed = average*(1 - beta) + max * beta;
	f_expansed = f(expansed.vec);
}

void DownhillSimplexMethod::Reduction() //�������� ����
{
	for (int i = 0; i < N + 1; i++)
	{
		smplx[i] = (smplx[i] - min)*0.5;
		F[i] = f(smplx[i].vec);
	}
}