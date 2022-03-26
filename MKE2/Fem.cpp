#include "Fem.h"

FEM::FEM()
{
   grid = new Grid();

   num_of_knots = grid->knots.size();
   num_of_FE = grid->fes.size();
   MakeSparseFormat();
   q = new real[num_of_knots]{};
   qprev = new real[num_of_knots]{};
   b = new real[num_of_knots]{};
   temp = new real[num_of_knots]{};
}

void FEM::MakeSparseFormat()
{
   const int N = 4;
   
   int* listbeg = new int[num_of_knots];//{-1};
   for (int i = 0; i < num_of_knots; i++)  listbeg[i] = -1;

    int* list1, * list2;
   list1 = new int[num_of_knots * num_of_knots]{};
   list2 = new int[num_of_knots * num_of_knots]{};

   int listsize = -1, iaddr, ind1, ind2, k;

   for (int iel = 0; iel < num_of_FE; iel++) // ï
   {
      for (int i = 0; i < N; i++) // 
      {
         k = grid->fes[iel]->knots_num[i]; //
         for (int j = i + 1; j < N; j++) // need to set N = ?
         {
            ind1 = k;
            ind2 = grid->fes[iel]->knots_num[j];  //
            if (ind2 < ind1) //
            {
               ind1 = ind2;
               ind2 = k;
            }
            iaddr = listbeg[ind2];
            if (iaddr == -1) // 
            {
               listsize++;
               listbeg[ind2] = listsize;
               list1[listsize] = ind1;
               list2[listsize] = -1;
            }
            else // 
            {
               while (list1[iaddr] < ind1 && list2[iaddr] >= 0)
                  iaddr = list2[iaddr];
               if (list1[iaddr] > ind1)  // 
               {                         // 
                  listsize++;
                  list1[listsize] = list1[iaddr];
                  list2[listsize] = list2[iaddr];
                  list1[iaddr] = ind1;
                  list2[iaddr] = listsize;
               }
               else if (list1[iaddr] < ind1) // 
               {
                  listsize++;
                  list2[iaddr] = listsize;
                  list1[listsize] = ind1;
                  list2[listsize] = -1;
               }
            }
         }
      }
   }

   A = new Matrix;
   A->ig = new int[num_of_knots + 1]{};
   A->jg = new int[listsize + 1]{};  // +1???

   for (int i = 0; i < num_of_knots; i++)
   {
      A->ig[i + 1] = A->ig[i];

      for (iaddr = listbeg[i]; iaddr != -1; )
      {
         A->jg[A->ig[i + 1]] = list1[iaddr];
         A->ig[i + 1]++;
         iaddr = list2[iaddr];
      }
   }

   delete[] listbeg;
   delete[] list1;
   delete[] list2;


   A->l = new real[listsize + 1]{};
   A->u = new real[listsize + 1]{};
   A->d = new real[num_of_knots]{};

}

void FEM::AddLocal(Matrix* A, long knot_num[4], real localA[4][4], real coeff)
{
   int ibeg, iend, ind;
   for (int i = 0; i < 4; i++)
      A->d[knot_num[i]] += coeff * localA[i][i];
   for (int i = 0; i < 4; i++)
   {
      ibeg = A->ig[knot_num[i]];
      for (int j = 0; j < i; j++)
      {
         iend = A->ig[knot_num[i] + 1]; 
         while (A->jg[ibeg] != knot_num[j])
         {
            ind = (ibeg + iend) / 2;
            if (A->jg[ind] < knot_num[j])
               ibeg = ind + 1;
            else
               iend = ind;
         }
         A->l[ibeg] += coeff * localA[i][j];
         A->u[ibeg] += coeff * localA[j][i];
         ibeg++;
      }
   }

}

void FEM::copy(real* from, real* to)
{
   for (int i = 0; i < num_of_knots; i++)
      to[i] = from[i];
}

void FEM::SolveLinear()
{
   CreateSLAE();
   SolveSLAE();

}

void FEM::Output(std::ofstream& out)
{
    out << endl << "Результат в узлах (веса):" << endl;
    out << " ___________________________________________________________________ " << endl;

    out.setf(ios::left);
    out.width(15);
    out << "| № элемента " << "  | ";
    out.width(15);
    out << "x" << "| ";
    out.width(15);
    out << "y" << "| ";
    out.width(15);
    out << "u*" << "|" << endl;
    out << "|----------------|----------------|----------------|----------------|" << endl;

    for (int i = 0; i < num_of_knots; i++)
    {
        out.unsetf(ios_base::floatfield);
        out << std::fixed;
        out.setf(ios::left);
        out << "| ";
        out.width(15);
        out << i + 1 << "| ";
        out.width(15);
        out << grid->knots[i]->x << "| ";
        out.width(15);
        out << grid->knots[i]->y << "| ";
        out.width(15);
        out.precision(6);
        out << scientific;
        out << q[i] << "| " << endl;
    }
}

void FEM::AddFirstBounds()
{
   for (auto cond : grid->bounds)
   {
      real C = 1e10; /// очень большое число
      A->d[cond->knot_num] = C;
      for (int j = A->ig[cond->knot_num]; j < A->ig[cond->knot_num + 1]; j++)
         A->l[j] = 0.;
      for (int j = 0; j < A->ig[num_of_knots]; j++)
         if (A->jg[j] == cond->knot_num)
            A->u[j] = 0.;

      b[cond->knot_num] = 0. * C;  
   }
}

void FEM::CreateSLAE()
{
   FE* fe;
   for (int i = 0; i < num_of_FE; i++)
   {
      fe = grid->fes[i];
      real localG[4][4];
      CreateLocalG(fe, localG);
      AddLocal(A, fe->knots_num, localG, 1);
      AddToB(fe);
   }

   for (int i = 0; i < num_of_knots; i++)
      std::cout << b[i] << '\n';
   //AddSecondBounds();
   AddFirstBounds();
}

void FEM::CreateLocalG(FE* fe, real localG[4][4])
{
   real lambda = 1./ grid->materials[fe->material - 1].mu; // mu - spline
   real mu = grid->materials[fe->material - 1].mu;
   real hx = grid->knots[fe->knots_num[1]]->x - grid->knots[fe->knots_num[0]]->x,
        hy = grid->knots[fe->knots_num[2]]->y - grid->knots[fe->knots_num[0]]->y;

   localG[0][0] = localG[1][1] = localG[2][2] = localG[3][3] =  lambda * (hy / hx + hx / hy) / 3.;
   localG[0][1] = localG[1][0] = localG[2][3] = localG[3][2] = -lambda * (hy / hx - hx / hy / 2.) / 3.;
   localG[0][2] = localG[2][0] = localG[1][3] = localG[3][1] = lambda * (hy / hx / 2. - hx / hy) / 3.;
   localG[0][3] = localG[3][0] = localG[1][2] = localG[2][1] = -lambda * (hy / hx + hx / hy) / 6.;
   //localG[0][0] = localG[1][1] = localG[2][2] = localG[3][3] =  (hy / hx + hx / hy) / 3. / mu;
   //localG[0][1] = localG[1][0] = localG[2][3] = localG[3][2] = -(hy / hx - hx / hy / 2.) / 3. / mu;
   //localG[0][2] = localG[2][0] = localG[1][3] = localG[3][1] =  (hy / hx / 2. - hx / hy) / 3. / mu;
   //localG[0][3] = localG[3][0] = localG[1][2] = localG[2][1] = -(hy / hx + hx / hy) / 6. / mu;
}

void FEM::AddToB(FE* fe)
{
   real f_[4]{};
   for (int i = 0; i < 4; i++)
      f_[i] = grid->materials[fe->material - 1].j;
   real hx = grid->knots[fe->knots_num[1]]->x - grid->knots[fe->knots_num[0]]->x,
        hy = grid->knots[fe->knots_num[2]]->y - grid->knots[fe->knots_num[0]]->y;

   real localb[4]{
   hx * hy * (4. * f_[0] + 2. * f_[1] + 2. * f_[2] +      f_[3]) / 36.,
   hx * hy * (2. * f_[0] + 4. * f_[1] +      f_[2] + 2. * f_[3]) / 36.,
   hx * hy * (2. * f_[0] +      f_[1] + 4. * f_[2] + 2. * f_[3]) / 36.,
   hx * hy * (     f_[0] + 2. * f_[1] + 2. * f_[2] + 4. * f_[3]) / 36.
   };

   for (int i = 0; i < 4; i++)
      b[fe->knots_num[i]] += localb[i];
}

real FEM::scalar(real* u, real* v, int size)
{
   real sum = 0;
   for (int i = 0; i < size; i++)
      sum += v[i] * u[i];

   return sum;
}

void FEM::MatxVec(real* v, Matrix* A, real* b)
{
   for (int i = 0; i < num_of_knots; i++)
      v[i] = A->d[i] * b[i];

   for (int i = 0; i < num_of_knots; i++)
      for (int j = A->ig[i]; j < A->ig[i + 1]; j++) // -1?
      {
         v[i] += A->l[j] * b[A->jg[j]];
         v[A->jg[j]] += A->u[j] * b[i];
      }
}

void FEM::SolveSLAE()
{
   real* z, * r, * p, * ff, * x;
   ff = new real[num_of_knots]{};
   z = new real[num_of_knots]{};
   r = new real[num_of_knots]{};
   p = new real[num_of_knots]{};
   real res, alpha, beta, skp, eps = 1e-17;
   int i, k;
   x = q;

   real lastres;
   MatxVec(ff, A, x);
   for (i = 0; i < num_of_knots; i++)
      z[i] = r[i] = b[i] - ff[i];
   MatxVec(p, A, z);
   res = sqrt(scalar(r, r, num_of_knots)) / sqrt(scalar(b, b, num_of_knots));

   for (k = 1; k < 100000 && res > eps; k++)
   {
      lastres = res;
      skp = scalar(p, p, num_of_knots);
      alpha = scalar(p, r, num_of_knots) / skp;
      for (i = 0; i < num_of_knots; i++)
      {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }
      MatxVec(ff, A, r);
      beta = -scalar(p, ff, num_of_knots) / skp;
      for (i = 0; i < num_of_knots; i++)
      {
         z[i] = r[i] + beta * z[i];
         p[i] = ff[i] + beta * p[i];
      }
      res = sqrt(scalar(r, r, num_of_knots)) / sqrt(scalar(b, b, num_of_knots));
   }
   //std::cout << "Residual: " << res << std::endl;
}