#include "EMproblem.h"

EMproblem::EMproblem()
{
	fem = new FEM();
	q = fem->q;

	GetB2 = [this](Knot point) {
		double x = point.x;
		double y = point.y;

		int igl[4]{};
		int i;

		Grid* grid = fem->grid;
		int num_fes = grid->fes.size();
		for (i = 0; i < num_fes; i++) // Определяем в каком конечном элементе точка
		{
			// Определяем глобальные узлы для локального элемента
			for (int j = 0; j < 4; j++)
				igl[j] = grid->fes[i]->knots_num[j];

			// Если указанные точки входят в область, то мы нашли, в каком элементе лежит точка
			if (x >= grid->knots[igl[0]]->x && x <= grid->knots[igl[1]]->x &&
				y >= grid->knots[igl[0]]->y && y <= grid->knots[igl[2]]->y)
				break;
		}

		double hx = grid->knots[igl[1]]->x - grid->knots[igl[0]]->x;
		double hy = grid->knots[igl[2]]->y - grid->knots[igl[0]]->y;

		real	X1 = X(0, point.x, grid->knots[igl[0]]->x, grid->knots[igl[1]]->x),
			X2 = X(1, point.x, grid->knots[igl[0]]->x, grid->knots[igl[1]]->x),
			Y1 = X(0, point.y, grid->knots[igl[0]]->y, grid->knots[igl[2]]->y),
			Y2 = X(1, point.y, grid->knots[igl[0]]->y, grid->knots[igl[2]]->y);

		double mes = hx * hy;

		double localG[4][4]{};
		fem->CreateLocalG(grid->fes[i], localG);

		double B2 = 0.;
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				B2 += grid->materials[grid->fes[i]->material - 1].mu * localG[j][k] * fem->q[igl[j]] * fem->q[igl[k]];
			}
		}
		B2 /= mes; // книжная версия

		//real b1, b2;
		//
		//b1 = (X1 * (fem->q[igl[2]] - fem->q[igl[0]]) + X2 * (fem->q[igl[3]] - fem->q[igl[1]])) / hy;
		//b2 = - (Y1 * (fem->q[igl[1]] - fem->q[igl[0]]) + Y2 * (fem->q[igl[3]] - fem->q[igl[2]])) / hx;
		//std::cout << b1 << " " << b2 << '\n';
		return B2;//b1 * b1 + b2 * b2;
	};

}


//std::function<real(Knot point)> EMproblem::GetB2(Knot point) // Получение значения в произвольной точке


real EMproblem::GetA(Knot point)
{
	double x = point.x;
	double y = point.y;

	int igl[4]{};
	int i;

	Grid* grid = fem->grid;
	int num_fes = grid->fes.size();
	for (i = 0; i < num_fes; i++) // Определяем в каком конечном элементе точка
	{
		// Определяем глобальные узлы для локального элемента
		for (int j = 0; j < 4; j++)
			igl[j] = grid->fes[i]->knots_num[j];

		// Если указанные точки входят в область, то мы нашли, в каком элементе лежит точка
		if (x >= grid->knots[igl[0]]->x && x <= grid->knots[igl[1]]->x &&
			y >= grid->knots[igl[0]]->y && y <= grid->knots[igl[2]]->y)
			break;
	}

	double hx = grid->knots[igl[1]]->x - grid->knots[igl[0]]->x;
	double hy = grid->knots[igl[2]]->y - grid->knots[igl[0]]->y;

	real	X1 = X(0, point.x, grid->knots[igl[0]]->x, grid->knots[igl[1]]->x),
		X2 = X(1, point.x, grid->knots[igl[0]]->x, grid->knots[igl[1]]->x),
		Y1 = X(0, point.y, grid->knots[igl[0]]->y, grid->knots[igl[2]]->y),
		Y2 = X(1, point.y, grid->knots[igl[0]]->y, grid->knots[igl[2]]->y);

	double u =
		fem->q[igl[0]] * X1 * Y1 +
		fem->q[igl[1]] * X2 * Y1 +
		fem->q[igl[2]] * X1 * Y2 +
		fem->q[igl[3]] * X2 * Y2;

	return u;
}

void EMproblem::Start()
{
	//fem->SolveLinear();
	fem->SolveNonlinear(GetB2);
	std::ofstream out("Result.txt");
	fem->Output(out);
	out.close();
}

void EMproblem::Output(std::ofstream& out)
{
	std::ifstream fpoints("Point");
	int n_point;
	fpoints >> n_point;
	out.setf(ios::left);
	out.width(15);
	out << "x";
	out.width(15);
	out << "y";
	out.width(15);
	out << "A";
	out.width(15);
	out << "|B|" << '\n';
	for (int i = 0; i < n_point; i++)
	{
		real x, y, A, B;
		fpoints >> x >> y;
		Knot point(x, y);
		A = GetA(point);
		B = GetB2(point);
		out.width(15);
		out << x;
		out.width(15);
		out << y;
		out.width(15);
		out << A; 
		out.width(15);
		out << sqrt(B) << '\n';
	}
	fpoints.close();
	
}
