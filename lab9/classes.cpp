#include "classes.h"

/*
References:
https://www.geeksforgeeks.org/area-of-a-polygon-with-given-n-ordered-vertices/

====================================================================================================
A big thanks to Milos Oundjian, he helped me a lot to debug this
*/

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator+=(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& a) {
    return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector& a, const Vector& b) {
    return Vector(a[0] / b[0], a[1] / b[1], a[2] / b[2]);
}
Vector operator/(const Vector& a, double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
bool operator<(const Vector& a, const Vector& b) {
    return (a[0] < b[0]) && (a[1] < b[1]) && (a[2] < b[2]);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector pow(const Vector& a, double b) {
    return Vector(std::min(std::pow(a[0], b),255.), std::min(std::pow(a[1], b), 255.), std::min(std::pow(a[2], b),255.));
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
double min(const Vector& a){
    return std::min(std::min(a[0], a[1]), a[2]);
}

double max(const Vector& a){
    return std::max(std::max(a[0], a[1]), a[2]);
}

double max_idx(const Vector& a){
    int idx = 0;
    for (size_t i = 0; i < 3; i++){
        if (a[i] > a[idx]){
            idx = i;
        }
    }
    return idx;
}

// Polygon member function definitions
// Debugging inspiration: https://github.com/ekaterina-borisova/graphics_cse306
double calculate_area_contribution(const Vector& current_vert, const Vector& next_vert) {
    return current_vert[0] * next_vert[1] - current_vert[1] * next_vert[0];
}

double sum_area_contributions(const std::vector<Vector>& vertices) {
    double poly_area = 0.0;
    size_t vertex_count = vertices.size();

    for (size_t idx = 0; idx < vertex_count; ++idx) {
        size_t next_idx = (idx < vertex_count - 1) ? (idx + 1) : 0;
        poly_area += calculate_area_contribution(vertices[idx], vertices[next_idx]);
    }

    return poly_area / 2.0;
}

double Polygon::compute_area() {
    size_t vertex_count = vertices.size();
    if (vertex_count < 3) return 0.0;

    double poly_area = sum_area_contributions(vertices);
    return std::abs(poly_area);
}


double calculate_triangle_area(const Vector& a, const Vector& b) {
    return std::abs(a[0] * b[1] - a[1] * b[0]) / 2.0;
}

double calculate_triangle_integral(const Vector triangle[3], const Vector& ref_point, double tri_area) {
    double integral_sum = 0.0;

    for (int k = 0; k < 3; ++k) {
        for (int l = k; l < 3; ++l) {
            integral_sum += tri_area / 6.0 * dot(triangle[k] - ref_point, triangle[l] - ref_point);
        }
    }

    return integral_sum;
}

double Polygon::distance_integral(Vector& ref_point) {
    double integral_sum = 0.0;
    size_t vertex_count = vertices.size();
    if (vertex_count < 3) return 0.0;

    for (size_t idx = 1; idx < vertex_count - 1; ++idx) {
        Vector triangle[3] = {vertices[0], vertices[idx], vertices[idx + 1]};
        Vector edge_a = triangle[1] - triangle[0];
        Vector edge_b = triangle[2] - triangle[0];
        double tri_area = calculate_triangle_area(edge_a, edge_b);

        integral_sum += calculate_triangle_integral(triangle, ref_point, tri_area);
    }

    return integral_sum;
}

Vector compute_centroid_contribution(const Vector& a, const Vector& b, double poly_area) {
    return (a + b) * (a[0] * b[1] - b[0] * a[1]) / (6.0 * poly_area);
}

Vector Polygon::calculate_centroid() {
    Vector centroid_vector(0.0, 0.0, 0.0);
    double poly_area = compute_area();
    size_t vertex_count = vertices.size();

    if (poly_area == 0.0) {
        return vertex_count > 0 ? vertices[0] : Vector(0.0, 0.0, 0.0);
    }

    for (size_t idx = 0; idx < vertex_count - 1; ++idx) {
        centroid_vector = centroid_vector + compute_centroid_contribution(vertices[idx], vertices[idx + 1], poly_area);
    }
    centroid_vector = centroid_vector + compute_centroid_contribution(vertices[vertex_count - 1], vertices[0], poly_area);

    return -centroid_vector;
}

