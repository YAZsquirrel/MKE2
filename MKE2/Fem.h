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
	int num_of_knots, num_of_FE, listsize;
	real* q, *qprev;
	Grid* grid;
	void CreateLocalG(FE* fe, real localG[4][4]);

private:

	struct Matrix
	{
		real* l, * u, * d;
		int* ig, * jg;
	};


	void MakeSparseFormat();  // �������

	// ��������� � ������� � ����� �� �����
	void AddLocal(Matrix* A, long knot_num[4], real localA[4][4], real coeff);


	void AddFirstBounds();

	void CreateSLAE();



	void AddToB(FE* fe);

	Matrix* A;
	Matrix* M;
	real* b, * temp;

	real f(Knot* knot_); //= Jz
	//real ug(Knot* knot_);

	void WriteMatrix(Matrix* A);
	void copy(real* x, real* y);
	real scalar(real* v, real* u, int size);
	void MatxVec(real* v, Matrix* A, real* b);
	void SolveSLAE();
	void CreateKholessky();
	void solve_L(real* f, real*& x);
	void solve_LT(real* f, real*& x);
	void solve_LLT(real* f, real*& x);

};
