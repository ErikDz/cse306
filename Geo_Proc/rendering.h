#ifndef RENDERING_H
#define RENDERING_H

#include <vector>
#include <string>
#include "utils/classes.h"

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0);

#endif // RENDERING_H
