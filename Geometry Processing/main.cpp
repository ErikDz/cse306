#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <cfloat>
#include <tuple>
#include <omp.h>
#include <list>
#include <chrono>
#include <iostream>

#include "svg.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./stb/stb_image_write.h"

// Initialize random number generator
static std::default_random_engine randomEngine(10);
static std::uniform_real_distribution<double> uniformDist(0.0, 1.0);

/* ------------------------ VECTOR AND POLYGON UTILITIES ------------------------ */

// Function to find the intersection point of a polygon edge with a vector
Vector findPolygonIntersection(const Vector& startPoint, const Vector& endPoint, const std::vector<Vector>& edge) {
    Vector normal = Vector(edge[1][1] - edge[0][1], edge[0][0] - edge[1][0], 0).normalized();
    double t = dot(normal, edge[0] - startPoint) / dot(normal, endPoint - startPoint);
    if (0 <= t && t <= 1) {
        return startPoint + t * (endPoint - startPoint);
    }
    return Vector();
}

// Function to check if a point is inside a polygon defined by its edges
bool isPointInsidePolygon(const Vector& point, const std::vector<Vector>& edge) {
    Vector normal = Vector(edge[1][1] - edge[0][1], edge[0][0] - edge[1][0], 0);
    return dot(normal, point - edge[0]) <= 0;
}

/* Polygon clipping algorithm using Sutherland-Hodgman algorithm */
Polygon clipSubjectPolygon(Polygon& subjectPolygon, Polygon& clippingPolygon) {
    Polygon outputPolygon;
    for (size_t i = 0; i < clippingPolygon.vertices.size(); ++i) {
        std::vector<Vector> clipEdge = {
            clippingPolygon.vertices[i],
            clippingPolygon.vertices[(i > 0) ? (i - 1) : (clippingPolygon.vertices.size() - 1)]
        };
        outputPolygon = Polygon();

        for (size_t j = 0; j < subjectPolygon.vertices.size(); ++j) {
            Vector currentVertex = subjectPolygon.vertices[j];
            Vector previousVertex = subjectPolygon.vertices[(j > 0) ? j - 1 : (subjectPolygon.vertices.size() - 1)];
            Vector intersection = findPolygonIntersection(previousVertex, currentVertex, clipEdge);

            if (isPointInsidePolygon(currentVertex, clipEdge)) {
                if (!isPointInsidePolygon(previousVertex, clipEdge)) {
                    outputPolygon.addVertex(intersection);
                }
                outputPolygon.addVertex(currentVertex);
            } else if (isPointInsidePolygon(previousVertex, clipEdge)) {
                outputPolygon.addVertex(intersection);
            }
        }
        subjectPolygon = outputPolygon;
    }
    return outputPolygon;
}

/* ------------------------ MAIN FUNCTION ------------------------ */

int main() {
    // Define a new subject polygon with its vertices
    Polygon subjectPolygon({
        Vector(0.2, 0.1), Vector(0.4, 0.7), Vector(0.6, 0.3),
        Vector(0.8, 0.8), Vector(0.7, 0.2), Vector(0.5, 0.5)
    });

    // Define a new clipping polygon with its vertices
    Polygon clippingPolygon({
        Vector(0.4, 0.4), Vector(0.4, 0.6),
        Vector(0.6, 0.6), Vector(0.6, 0.4)
    });

    // Save the initial polygons to an SVG file
    save_svg({subjectPolygon, clippingPolygon}, "images/new_initial.svg", "none");

    // Perform the polygon clipping and save the result to an SVG file
    save_svg({clipSubjectPolygon(subjectPolygon, clippingPolygon)}, "images/new_clipped.svg", "none");

    return 0;
}
