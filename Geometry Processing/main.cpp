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

/*------------------------- Vornoi + vornoi utils ----------------------------------------------*/

// Function to calculate the normal vector for Voronoi diagram
Vector calculateNormal(const Vector& pointA, const Vector& pointB, double weightA, double weightB) {
    return (pointA + pointB) / 2 + (weightA - weightB) * (pointB - pointA) / (2 * std::pow((pointA - pointB).squaredNorm(), 2));
}

// Function to find the intersection point for the Voronoi power diagram
Vector voronoiIntersection(
    const Vector& pointA,
    const Vector& pointB,
    const std::vector<Vector>& edgePoints,
    const std::vector<double>& edgeWeights
) {
    Vector normal = calculateNormal(edgePoints[0], edgePoints[1], edgeWeights[0], edgeWeights[1]);
    double t = dot(normal - pointA, edgePoints[0] - edgePoints[1]) / dot(pointB - pointA, edgePoints[0] - edgePoints[1]);
    if (0 <= t && t <= 1) {
        return pointA + t * (pointB - pointA);
    }
    return Vector();
}

// Function to check if a point is inside the Voronoi cell
bool isPointInsideVoronoi(const Vector& point, const std::vector<Vector>& edgePoints, const std::vector<double>& edgeWeights) {
    Vector normal = calculateNormal(edgePoints[0], edgePoints[1], edgeWeights[0], edgeWeights[1]);
    return dot(point - normal, edgePoints[1] - edgePoints[0]) < 0;
}

/* Voronoi Parallel Linear Enumeration */
std::vector<Polygon> voronoiParallelLinearEnumeration(
    const std::vector<Vector>& sites,
    const Polygon& boundingPolygon,
    const std::vector<double>& siteWeights
) {
    std::vector<Polygon> voronoiCells(sites.size());

    #pragma omp parallel for
    for (int i = 0; i < sites.size(); i++) {
        Vector currentSite = sites[i];
        Polygon currentEdges = boundingPolygon;

        for (int j = 0; j < sites.size(); j++) {
            if (i != j) {
                Vector otherSite = sites[j];
                std::vector<Vector> edgePoints = {currentSite, otherSite};
                std::vector<double> edgeWeights = {siteWeights[i], siteWeights[j]};
                Polygon clippedPolygon;

                for (int k = 0; k < currentEdges.vertices.size(); k++) {
                    Vector currentVertex = currentEdges.vertices[k];
                    Vector previousVertex = currentEdges.vertices[(k > 0) ? (k - 1) : (currentEdges.vertices.size() - 1)];
                    Vector intersection = voronoiIntersection(previousVertex, currentVertex, edgePoints, edgeWeights);

                    if (isPointInsideVoronoi(currentVertex, edgePoints, edgeWeights)) {
                        if (!isPointInsideVoronoi(previousVertex, edgePoints, edgeWeights)) {
                            clippedPolygon.addVertex(intersection);
                        }
                        clippedPolygon.addVertex(currentVertex);
                    } else if (isPointInsideVoronoi(previousVertex, edgePoints, edgeWeights)) {
                        clippedPolygon.addVertex(intersection);
                    }
                }
                currentEdges = clippedPolygon;
            }
        }
        voronoiCells[i] = currentEdges;
    }
    return voronoiCells;
}

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
    if (false){
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
    }

    int n_points = 10;
    Polygon boundingPolygon({
       Vector(0.,0.) , Vector(0., 1.), Vector(1., 1.), Vector(1., 0.)
    });
    std::vector<Vector> sites(n_points);
    for (int i=0; i<n_points; i++) {
        sites[i] = Vector(uniformDist(randomEngine), uniformDist(randomEngine));
    }
    save_svg(voronoiParallelLinearEnumeration(sites, boundingPolygon, std::vector<double>(n_points, 1.)), "images/voronoi.svg", "none");

    return 0;
}
