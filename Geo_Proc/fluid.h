#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include "utils/classes.h"
#include "./liblbfgs/lbfgs.h"
#include "voronoi.h"  


std::vector<Polygon> construct_fluid(const Vector* points, const double* weights, int num_points);
extern "C" lbfgsfloatval_t eval_f(void *instance, const lbfgsfloatval_t *weights, lbfgsfloatval_t *gradient, const int num_points, const lbfgsfloatval_t step);

#endif // FLUID_H
