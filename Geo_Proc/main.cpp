#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <thread>
#include <numeric>
#include <iterator>
#include <random>

#include "svg.h"
#include "fluid.h"
#include "utils/classes.h"
#include "voronoi.h"
#include "fluid.h"
#include "rendering.h"
#include "gallouet.h"
#include "./liblbfgs/lbfgs.c"

static std::default_random_engine engine(1);
static std::uniform_real_distribution<double> uniform (0 ,1);



// To generate the Voronoi diagrams:
/*
int main() {
    int num_points = 100;
    Vector* points = new Vector[num_points];
    double* weights = new double[num_points];

    // Generate random points and weights
    for (int i = 0; i < num_points; i++) {
        points[i] = Vector(uniform(engine), uniform(engine));
        weights[i] = 0;
    }

    // Compute the Voronoi diagram
    std::vector<Polygon> diagram = generate_diagram(points, weights, num_points);

    // Save the Voronoi diagram as an SVG file
    save_svg(diagram, "images/voronoi_diagram.svg", "none");

    // Generate random points and weights
    for (int i = 0; i < num_points; i++) {
        points[i] = Vector(uniform(engine), uniform(engine));
        weights[i] = uniform(engine); 
    }

    // Compute the power diagram
    std::vector<Polygon> power_diagram = generate_diagram(points, weights, num_points);

    // Save the power diagram as an SVG file
    save_svg(power_diagram, "images/power_diagram.svg", "none");

    // Clean up dynamically allocated memory
    delete[] points;
    delete[] weights;


    return 0;
}*/
//For fluid simulation:
int main() {
    int point_count = 20;
    std::vector<Vector> positions(point_count);
    std::vector<Vector> velocities(point_count);
    std::vector<double> weights(point_count, 0);

    for (int i = 0; i < point_count; ++i) {
        positions[i][0] = rand() / static_cast<double>(RAND_MAX);
        velocities[i][0] = 0.0;
        positions[i][1] = rand() / static_cast<double>(RAND_MAX);
        velocities[i][1] = 0.0;

        velocities[i][2] = 0.0;
        positions[i][2] = 0.0;
        weights[i] = 0.05;
    }

    for (int time = 0; time < 100; ++time) {
        gal_step(positions, velocities, weights, time);
    }

    return 0;
}