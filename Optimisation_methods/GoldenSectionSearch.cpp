#include "GoldenSectionSearch.h"


GoldenSectionSearch::GoldenSectionSearch()
{
}


GoldenSectionSearch::~GoldenSectionSearch()
{
}

pair<double, double> GoldenSectionSearch::FindInterval(double lambda0, double d) 
{
	double lambda1, lambda2;
	if (f((x + S*lambda0).vec) < f((x + S*(lambda0 + d)).vec))
		d = -d;
	do
	{
		lambda1 = lambda0 + d;
		d *= 2;
		lambda2 = lambda1 + d;
		lambda0 = lambda1;

	} while (f((x + S*lambda1).vec) > f((x + S*(lambda2)).vec));
	lambda0 -= d / 2;
	return make_pair(min(lambda0, lambda2), max(lambda0, lambda2));
}

double GoldenSectionSearch::GSS(func _f, Vertex &_x, Vertex &_S)
{
	S = _S;	x = _x;	f = _f;
	auto interval = FindInterval(1, 0.1);
	auto l1 = interval.first + 0.381966011*(interval.second - interval.first);
	auto l2 = interval.second - 0.381966011*(interval.second - interval.first);
	auto f_l1 = f((x + S*l1).vec);
	auto f_l2 = f((x + S*l2).vec);
	while (abs(interval.second - interval.first) > 1e-13)
		if (f_l1 < f_l2)
		{
			interval.second = l2;
			l2 = l1;
			f_l2 = f_l1;
			l1 = interval.first + 0.381966011*(interval.second - interval.first);
			f_l1 = f((x + S*l1).vec);
		}
		else
		{
			interval.first = l1;
			l1 = l2;
			f_l1 = f_l2;
			l2 = interval.second - 0.381966011*(interval.second - interval.first);
			f_l2 = f((x + S*l2).vec);
		}
	return (interval.first + interval.second) / 2;
}
