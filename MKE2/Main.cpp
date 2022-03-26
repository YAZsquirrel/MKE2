#include "EMproblem.h"
#include "locale.h"
int main()
{
	setlocale(LC_ALL,"RU");

	EMproblem* emp = new EMproblem();
	emp->Start();
	std::ofstream out("ResAB.txt");
	emp->Output(out);


	return 0;
}