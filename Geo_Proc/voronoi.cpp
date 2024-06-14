#include "voronoi.h"
#include <cmath>

/*
References:
https://courses.cs.washington.edu/courses/cse326/00wi/projects/voronoi.html
Took some code from: https://www.geeksforgeeks.org/voronoi-diagram/ and https://github.com/JCash/voronoi
For debugging browsed:
https://github.com/jbhlevy/CSE306/blob/main/Project_2/main.cpp
*/

Vector calculate_normal_vector(const Vector& start, const Vector& end) {
    return Vector(end[1] - start[1], start[0] - end[0]);
}

Vector calculate_intersect_point(const Vector& prev_vert, const Vector& curr_vert, const Vector& start, const Vector& normal_vec) {
    double t_val = dot(start - prev_vert, normal_vec) / dot(curr_vert - prev_vert, normal_vec);
    return prev_vert + t_val * (curr_vert - prev_vert);
}

bool is_inside(const Vector& point, const Vector& start, const Vector& normal_vec) {
    return dot(start - point, normal_vec) <= 0;
}

Polygon edge_clipping(Polygon& poly, const Vector& start, const Vector& end) {
    Polygon newPoly;
    int vert_count = poly.vertices.size();
    Vector normal_vec = calculate_normal_vector(start, end);

    for (int idx = 0; idx < vert_count; idx++) {
        const Vector& curr_vert = poly.vertices[idx];
        const Vector& prev_vert = poly.vertices[(idx > 0) ? (idx - 1) : (vert_count - 1)];
        Vector intersect_pt = calculate_intersect_point(prev_vert, curr_vert, start, normal_vec);

        if (is_inside(curr_vert, start, normal_vec)) {
            if (!is_inside(prev_vert, start, normal_vec)) {
                newPoly.vertices.push_back(intersect_pt);
            }
            newPoly.vertices.push_back(curr_vert);
        } else if (is_inside(prev_vert, start, normal_vec)) {
            newPoly.vertices.push_back(intersect_pt);
        }
    }

    return newPoly;
}

Vector calculate_mid_point(const Vector& pt1, const Vector& pt2, double wt1, double wt2) {
    Vector mid_pt = (pt1 + pt2) / 2;
    Vector delta = pt2 - pt1;
    mid_pt = mid_pt + (wt1 - wt2) / (2 * (-delta).squaredNorm()) * delta;
    return mid_pt;
}

Vector calculate_voronoi_intersect(const Vector& prev_vert, const Vector& curr_vert, const Vector& mid_pt, const Vector& delta) {
    double t_val = dot(mid_pt - prev_vert, delta) / dot(curr_vert - prev_vert, delta);
    return prev_vert + t_val * (curr_vert - prev_vert);
}

bool is_closer_to_point(const Vector& point, const Vector& pt1, const Vector& pt2, double wt1, double wt2) {
    return (point - pt1).squaredNorm() - wt1 <= (point - pt2).squaredNorm() - wt2;
}

Polygon clip_voronoi(Polygon& poly, int idx1, int idx2, const Vector* pts, const double* wts) {
    Polygon clipped_poly;
    int vert_count = poly.vertices.size();

    for (int vert_idx = 0; vert_idx < vert_count; vert_idx++) {
        const Vector& curr_vert = poly.vertices[vert_idx];
        const Vector& prev_vert = poly.vertices[(vert_idx > 0) ? (vert_idx - 1) : (vert_count - 1)];
        Vector mid_pt = calculate_mid_point(pts[idx1], pts[idx2], wts[idx1], wts[idx2]);
        Vector delta = pts[idx2] - pts[idx1];
        Vector intersect_pt = calculate_voronoi_intersect(prev_vert, curr_vert, mid_pt, delta);

        if (is_closer_to_point(curr_vert, pts[idx1], pts[idx2], wts[idx1], wts[idx2])) {
            if (!is_closer_to_point(prev_vert, pts[idx1], pts[idx2], wts[idx1], wts[idx2])) {
                clipped_poly.vertices.push_back(intersect_pt);
            }
            clipped_poly.vertices.push_back(curr_vert);
        } else if (is_closer_to_point(prev_vert, pts[idx1], pts[idx2], wts[idx1], wts[idx2])) {
            clipped_poly.vertices.push_back(intersect_pt);
        }
    }

    return clipped_poly;
}

// https://github.com/steven-mi/voronoi-diagram/blob/master/src/generator/VoronoiGenerator.cpp
std::vector<Polygon> generate_diagram(const Vector* pt_arr, const double* wt_arr, int pt_count) {
    std::vector<Polygon> poly_diagram(pt_count);
    
    for (int idx = 0; idx < pt_count; idx++) {
        poly_diagram[idx].vertices = {Vector(0, 0, 0), Vector(0, 1, 0), Vector(1, 1, 0), Vector(1, 0, 0)};

        for (int other_idx = 0; other_idx < pt_count; other_idx++) {
            if (idx != other_idx) {
                poly_diagram[idx] = clip_voronoi(poly_diagram[idx], idx, other_idx, pt_arr, wt_arr);
            }
        }
    }

    return poly_diagram;
}
