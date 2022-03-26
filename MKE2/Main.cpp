#include "Fem.h"
#include "locale.h"
int main()
{
	setlocale(LC_ALL,"RU");
	FEM *fem = new FEM();
	fem->SolveLinear();
	std::ofstream out("Result.txt");
	fem->Output(out);

	return 0;
}