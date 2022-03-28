#pragma once
#include "Grid.h"
#include "Fem.h"
#include <functional>
//using namespace FEM_ns;


class EMproblem
{
	public:
	EMproblem();
	FEM *fem;

	real* q;
	std::function<real(Knot point)> GetB2;
	real GetA(Knot point);
	void Start();
	void Output(std::ofstream& out);
	private: 
	std::function<real(int i, real x, real x1, real x2)> X = [](int i, real x, real x1, real x2)
	{
		if (i % 2 == 0) 
			return (x2 - x) / (x2 - x1);
		else 
			return (x - x1) / (x2 - x1);
	};
};

