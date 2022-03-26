#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "math.h"
typedef double real;
using namespace std;

struct Knot
   {
      Knot(real _x, real _y) : x(_x), y(_y) {}
      Knot() : x(0.0), y(0.0) {}
      real x, y;
   };

   struct FE 
   {
      long knots_num[4]{};
      long material;
   };

   struct Bound 
   {
      long knot_num;
   };

   struct Material
   {
      real j;
      real mu;
   };

class Grid
{
public:
    vector<Knot*> knots;        // Координаты узлов
    vector<double> x;           // Сетка по x
    vector<double> y;           // Сетка по y
    vector<FE*> fes;        // Конечные элементы
    vector<Bound*> bounds;      // Границы
    vector<Material> materials; // Материалы
    Grid();

private:
    int CreateXY();
    int ReadCountElements(int& countKnots, int& countFEs, int& countBounds, string pathFile);
    int ReadCoordinates(string pathFile);
    int ReadKEs(string pathFile);
    int ReadBounds(string pathFile);
    int ReadMaterials(string pathFileJ, string pathFileMu);
    int ReadFEs(string pathFileMaterial, string pathFileNumsKnots);
};
