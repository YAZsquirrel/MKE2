#include "EMproblem.h"

//real* EMproblem::GetBInPoint(Knot point, FE* fe, real* q)
//{  
//    real hx = knots[fe->knots_num[1]].x - knots[fe->knots_num[0]].x,
//         hy = knots[fe->knots_num[2]].x - knots[fe->knots_num[0]].x;
//
//    real X1 = X(0, point.x, knots[fe->knots_num[0]].x, knots[fe->knots_num[1]].x),
//         X2 = X(1, point.x, knots[fe->knots_num[0]].x, knots[fe->knots_num[1]].x),
//         Y1 = X(0, point.y, knots[fe->knots_num[0]].y, knots[fe->knots_num[2]].y),
//         Y2 = X(1, point.y, knots[fe->knots_num[0]].y, knots[fe->knots_num[2]].y);
//
//    real* B = new real[2]{(X1 * (q[fe->knots_num[2]] - q[fe->knots_num[0]]) +
//                           X2 * (q[fe->knots_num[3]] - q[fe->knots_num[1]])) / hy ,
//                          (Y1 * (q[fe->knots_num[1]] - q[fe->knots_num[0]]) +
//                           Y2 * (q[fe->knots_num[3]] - q[fe->knots_num[2]])) / hx };
//    return B;
//}
//
//real EMproblem::GetAInPoint(Knot point, FE* fe, real* q)
//{
//   real X1 = X(0, point.x, knots[fe->knots_num[0]].x, knots[fe->knots_num[1]].x),
//        X2 = X(1, point.x, knots[fe->knots_num[0]].x, knots[fe->knots_num[1]].x),
//        Y1 = X(0, point.y, knots[fe->knots_num[0]].y, knots[fe->knots_num[2]].y),
//        Y2 = X(1, point.y, knots[fe->knots_num[0]].y, knots[fe->knots_num[2]].y);
//   return X1 * Y1 * q[fe->knots_num[0]] + X2 * Y1 * q[fe->knots_num[1]] + X1 * Y2 * q[fe->knots_num[2]] + X2 * Y2 * q[fe->knots_num[3]];
//}
