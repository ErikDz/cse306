#ifndef VORONOI_H
#define VORONOI_H

#include <vector>
#include "utils/classes.h"

Polygon edge_clipping(Polygon& poly, const Vector& start, const Vector& end);
Polygon clip_voronoi(Polygon& poly, int idx1, int idx2, const Vector* pts, const double* wts); 
std::vector<Polygon> generate_diagram(const Vector* pt_arr, const double* wt_arr, int pt_count);

#endif // VORONOI_H
