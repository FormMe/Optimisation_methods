#pragma once
#include "QuasiNewtonMethod.h"

class GoldfarbMethod : public QuasiNewtonMethod
{
public:
	GoldfarbMethod(ifstream &fin);
	~GoldfarbMethod();
private:
	virtual void ApproachMatrix() override;
};

