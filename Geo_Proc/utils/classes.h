#ifndef CLASSES_H
#define CLASSES_H

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>

// Vector class representing a 3D vector
class Vector
{
public:
    explicit Vector(double x = 0, double y = 0, double z = 0)
    {
        components[0] = x;
        components[1] = y;
        components[2] = z;
    }
    double squaredNorm() const
    {
        return components[0] * components[0] + components[1] * components[1] + components[2] * components[2];
    }
    double norm() const
    {
        return sqrt(squaredNorm());
    }
    Vector normalized()
    {
        double n = norm();
        Vector normalizedVector;
        normalizedVector[0] = components[0] / n;
        normalizedVector[1] = components[1] / n;
        normalizedVector[2] = components[2] / n;
        return normalizedVector;
    }
    double operator[](int i) const { return components[i]; }
    double &operator[](int i) { return components[i]; }
    double components[3];
};

// Function declarations
Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator-(const Vector& a);
Vector operator*(double a, const Vector& b);
Vector operator*(const Vector& a, double b);
Vector operator*(const Vector& a, const Vector& b);
Vector operator/(const Vector& a, const Vector& b);
Vector operator/(const Vector& a, double b);
bool operator<(const Vector& a, const Vector& b);
double dot(const Vector& a, const Vector& b);
Vector pow(const Vector& a, double b);
Vector cross(const Vector& a, const Vector& b);
double min(const Vector& a);
double max(const Vector& a);
double max_idx(const Vector& a);

class Polygon {
public:
    std::vector<Vector> vertices;

    double compute_area();
    double distance_integral(Vector& point);
    Vector calculate_centroid();
};

#endif // CLASSES_H
