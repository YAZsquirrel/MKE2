#include "Grid.h"
#define PI 3.14159265358979

//#define _CRT_SECURE_NO_WARNINGS
Grid::Grid()
{
	materials.resize(4);
	if (ReadMaterials("toku", "mu")) exit(1);

	int countKnots, countFEs, countBounds;
	if (ReadCountElements(countKnots, countFEs, countBounds, "inf2tr.dat")) exit(2);
	fes.resize(countFEs);
	knots.resize(countKnots);
	bounds.resize(countBounds);

	if (ReadCoordinates("rz.dat")) exit(3);
	CreateXY();
	if (ReadBounds("l1.dat")) exit(4);
	if (ReadFEs("nvkat2d.dat", "nvtr.dat")) exit(5);
}

int Grid::CreateXY()
{
	x.push_back(knots[0]->x);


	int countKnots = knots.size();
	double xLast = knots[1]->x;
	int i;
	for (i = 1; i < countKnots && xLast != x[0]; i++, xLast = knots[i]->x)		x.push_back(xLast);

	for (int j = 1; j < countKnots; j+=i)	
		y.push_back(knots[j]->y);

	return 0;
}

int Grid::ReadCountElements(int& countKnots, int& countFEs, int& countBounds, string pathFileCounts)
{
	ifstream in(pathFileCounts);

	if (in.is_open())
	{
		string line;
		getline(in, line);	// ������ ������ ������ �����

		if (!in.eof())
		{
			in >> line >> countKnots >> line >> countFEs >> line >> countBounds;
		}
		else
		{
			cout << "������ ������ ���������� ���������";
			return 2;
		}
	}
	else
	{
		cout << "������ �������� �����: " << pathFileCounts;
		return 1;
	}

	return 0;
}

int Grid::ReadCoordinates(string pathFileCoords)
{
	FILE* binInCoords; 
	fopen_s(&binInCoords, pathFileCoords.c_str(), "rb");
	if (!binInCoords)
	{
		std::cout << "������ �������� �����: " << pathFileCoords;
		return 1;
	}

	int nReadGood = 0; 
	int countKnots = knots.size();
	int r;
	for (int i = 0; i < countKnots; i++)
	{
		Knot* knot = knots[i] = new Knot;
		r = fread(knot, sizeof(*knot), 1, binInCoords); if (r < 1) break;  
		nReadGood += r * sizeof(*knot) / sizeof(double);
	}

	fclose(binInCoords);
	if (nReadGood / 2 != countKnots)
	{
		std::cout << "������ ������: ����������� ���������� ��������� �� ������������� ���������� ������";
		return(2);
	}

	return 0;
}

int Grid::ReadBounds(string pathFileBounds)
{
	FILE* binInBounds;
	fopen_s(&binInBounds, pathFileBounds.c_str(), "rb");


	if (!binInBounds)
	{
		std::cout << "������ �������� �����: " << pathFileBounds;
		return 1;
	}

	int nReadGood = 0; 
	int countBounds = bounds.size();
	int r;
	for (int i = 0; i < countBounds; i++)
	{
		Bound* bound = bounds[i] = new Bound;
		r = fread(bound, sizeof(*bound), 1, binInBounds);  if (r < 1 || bound->knot_num < 0) break;
		nReadGood++;
		bound->knot_num--;
	}

	fclose(binInBounds);
	if (nReadGood != countBounds)
	{
		std::cout << "������ ������: ����������� ���������� ������ �� ������������� ���������� ������ � ������� ��������";
		return(2);
	}

	return 0;

}

int Grid::ReadMaterials(string pathFileJ, string pathFileMu)
{
	ifstream inJ(pathFileJ);
	ifstream inMu(pathFileMu);

	if (inJ.is_open() || inMu.is_open())
	{
		long NumMaterial = 1;
		while (!inJ.eof() || !inMu.eof())
		{
			inJ >> NumMaterial;
			inJ >> materials[NumMaterial - 1].j;
			inMu >> NumMaterial;
			inMu >> materials[NumMaterial-1].mu;
			materials[NumMaterial - 1].mu *= 4.0 * PI * 10e-7;
		}
	}
	else
	{
		std::cout << "������ �������� ������ ��� ������ ����������.";
		return 1;
	}

	return 0;
}

int Grid::ReadFEs(string pathFileMaterial, string pathFileNumsKnots)
{
	FILE* binInMaterial;
	fopen_s(&binInMaterial, pathFileMaterial.c_str(), "rb");

	if (!binInMaterial)
	{
		std::cout << "������ �������� �����: " << pathFileMaterial;
		return 1;
	}

	FILE* binInNumsKnots;
	fopen_s(&binInNumsKnots, pathFileNumsKnots.c_str(), "rb");
	if (!binInNumsKnots)
	{
		std::cout << "������ �������� �����: " << pathFileNumsKnots;
		return 1;
	}

	struct recordNumsKnots { int num0, num1, num2, num3, ign0, ign1; } rec;

	int nReadGood = 0; // ����� ������� ����������� �����
	int countFE = fes.size();
	int r;
	for (int i = 0; i < countFE; i++)
	{
		FE* fe = fes[i] = new FE;

		r = fread(&(fe->material), sizeof(fe->material), 1, binInMaterial);  if (r < 1 || fe->material <= 0) break;
		r = fread(&rec, sizeof(rec), 1, binInNumsKnots); if (r < 1) break;  // ��������, ��� ������ ������ �������
		if (rec.num0 < 0 || rec.num1 < 0|| rec.num2 < 0 ||rec.num3 < 0) break;
		fe->knots_num[0] = rec.num2-1;
		fe->knots_num[1] = rec.num3-1;
		fe->knots_num[2] = rec.num0-1;
		fe->knots_num[3] = rec.num1-1;

		nReadGood++;
	}

	fclose(binInMaterial);
	fclose(binInNumsKnots);
	if (nReadGood != countFE)
	{
		std::cout << "������ ������: ����������� ���������� ��������� �� ������������� ���������� �������� ���������";
		return(2);
	}

	return 0;
}

