#ifndef GALLOUET_H
#define GALLOUET_H

#include <vector>
#include "utils/classes.h"

void gal_step(std::vector<Vector> &points, std::vector<Vector> &velocities, std::vector<double> &weights, int time_step);

#endif // GALLOUET_H
