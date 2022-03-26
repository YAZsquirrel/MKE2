#pragma once
#include "Grid.h"
#include "Fem.h"
#include <functional>
//using namespace FEM_ns;

class EMproblem
{
	public:
	FEM *fem = new FEM();

	real* q;
	real* GetBInPoint(Knot point, FE* fe, real* q);
	real GetAInPoint(Knot point, FE* fe, real* q);

	private: 
	std::function<real(int i, real x, real x1, real x2)> X = [](int i, real x, real x1, real x2)
	{
		if (i % 2 == 0) 
			return (x2 - x) / (x2 - x1);
		else 
			return (x - x1) / (x2 - x1);
	};
};

