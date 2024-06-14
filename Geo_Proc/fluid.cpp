#include "fluid.h"
#include <cmath>

// Define the compute_hessian function

/*
References:

Read when debugging:
https://github.com/shayyn20/CSE306/blob/master/assignment%202/func/fluidSimulation.h
*/
Polygon create_fluid_circle(const Vector &center, double radius_val, int segment_count) {
    Polygon fluid_polygon;
    fluid_polygon.vertices.resize(segment_count);
    std::vector<Vector> vertices(segment_count);

    for (int seg_idx = 0; seg_idx < segment_count; seg_idx++) {
        double angle_val = static_cast<double>(seg_idx) / segment_count * 2.0 * M_PI;
        vertices[seg_idx] = Vector(std::cos(angle_val), -std::sin(angle_val)) * radius_val + center;
    }

    fluid_polygon.vertices = vertices;
    return fluid_polygon;
}


void clip_polygon_with_points(Polygon &polygon, int idx, const Vector* point_set, const double* weight_set, int point_count) {
    for (int clip_idx = 0; clip_idx < point_count; clip_idx++) {
        if (idx != clip_idx) {
            polygon = clip_voronoi(polygon, idx, clip_idx, point_set, weight_set);
        }
    }
}

void clip_polygon_edges(Polygon &polygon, const std::vector<Vector> &vertices, int segment_count) {
    for (int vert_idx = 0; vert_idx < segment_count - 1; vert_idx++) {
        polygon = edge_clipping(polygon, vertices[vert_idx], vertices[vert_idx + 1]);
    }
    polygon = edge_clipping(polygon, vertices[segment_count - 1], vertices[0]);
}

std::vector<Polygon> construct_fluid(const Vector* point_set, const double* weight_set, int point_count) {
    std::vector<Polygon> fluid_diagram_set(point_count - 1);
    constexpr int segment_count = 200;

    #pragma omp parallel for 
    for (int idx = 0; idx < point_count - 1; idx++) {
        double radius_val = std::sqrt(weight_set[idx] - weight_set[point_count - 1]);
        Polygon fluid_polygon = create_fluid_circle(point_set[idx], radius_val, segment_count);

        fluid_diagram_set[idx].vertices = {Vector(0, 0, 0), Vector(0, 1, 0), Vector(1, 1, 0), Vector(1, 0, 0)};

        clip_polygon_with_points(fluid_diagram_set[idx], idx, point_set, weight_set, point_count);
        clip_polygon_edges(fluid_diagram_set[idx], fluid_polygon.vertices, segment_count);
    }

    return fluid_diagram_set;
}


lbfgsfloatval_t calculate_total_cost(
    std::vector<Polygon> &fluid_diagram_set,
    const lbfgsfloatval_t *weight_set,
    lbfgsfloatval_t *grad_set,
    Vector *point_set,
    const int point_count,
    double desired_area,
    double &total_fluid_area_sum
) {
    lbfgsfloatval_t total_cost = 0.0;

    for (int idx = 0; idx < point_count - 1; idx++) {
        double polygon_area_val = fluid_diagram_set[idx].compute_area();
        total_fluid_area_sum += polygon_area_val;
        grad_set[idx] = polygon_area_val - desired_area;
        total_cost += -fluid_diagram_set[idx].distance_integral(point_set[idx]) 
                      - weight_set[idx] * desired_area 
                      + weight_set[idx] * polygon_area_val;
    }

    return total_cost;
}

lbfgsfloatval_t calculate_final_cost(
    lbfgsfloatval_t total_cost,
    const lbfgsfloatval_t *weight_set,
    lbfgsfloatval_t *grad_set,
    const int point_count,
    double total_fluid_area_sum,
    double air_fraction_const
) {
    grad_set[point_count - 1] = 1.0 - total_fluid_area_sum - air_fraction_const;
    total_cost += weight_set[point_count - 1] * (1.0 - total_fluid_area_sum) - weight_set[point_count - 1] * air_fraction_const;

    return total_cost;
}

// Used extern "C" https://learn.microsoft.com/fr-fr/cpp/cpp/extern-cpp?view=msvc-170 
// Due to problems with import
// Got the suggestion to fix my bug from a stackoverflow post (lost it), this one helped implement: https://stackoverflow.com/questions/20737987/extern-c-when-exactly-to-use/20738072#20738072
// I apologise if this is not the best way to fix the problem; couldn't find a better solution
/**/
extern "C" lbfgsfloatval_t eval_f(
    void *point_instance,
    const lbfgsfloatval_t *weight_set,
    lbfgsfloatval_t *grad_set,
    const int point_count,
    const lbfgsfloatval_t step_size
) {
    Vector* point_set = static_cast<Vector*>(point_instance);
    std::vector<Polygon> fluid_diagram_set = construct_fluid(point_set, weight_set, point_count);
    
    double total_fluid_area_sum = 0.0;
    const double fluid_fraction_const = 0.4; //0.4
    const double air_fraction_const = 0.6; //0.6
    const double desired_area = fluid_fraction_const / (point_count-1);
    lbfgsfloatval_t total_cost = calculate_total_cost(
        fluid_diagram_set, weight_set, grad_set, point_set, point_count, desired_area, total_fluid_area_sum
    );

    total_cost = calculate_final_cost(total_cost, weight_set, grad_set, point_count, total_fluid_area_sum, air_fraction_const);

    return total_cost;
}
