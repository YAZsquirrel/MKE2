#pragma once
#include <list>
#include "Grid.h"

typedef double real;


class FEM
{

public:
	FEM();
	void SolveNonlenear();
	void SolveLinear();
	void Output(std::ofstream& out);

private:

	struct Matrix
	{
		real* l, * u, * d;
		int* ig, * jg;
	};

	Grid* grid;

	void MakeSparseFormat();  // Портрет

	// Добавляет в матрицу А какие то штуки
	void AddLocal(Matrix* A, long knot_num[4], real localA[4][4], real coeff);


	void AddFirstBounds();

	void CreateSLAE();

	void CreateLocalG(FE* fe, real localG[4][4]);


	void AddToB(FE* fe);

	Matrix* A;
	int num_of_knots, num_of_FE;
	real* b, * temp;
	real* q, *qprev;

	real f(Knot* knot_); //= Jz
	//real ug(Knot* knot_);

	void copy(real* x, real* y);
	real scalar(real* v, real* u, int size);
	void MatxVec(real* v, Matrix* A, real* b);
	void SolveSLAE();


};
