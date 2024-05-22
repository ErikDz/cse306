#include <cmath>
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
    double operator[](int i) const { return components[i]; };
    double &operator[](int i) { return components[i]; };
    double components[3];
};

// Vector arithmetic operations, I credit my Linear Algebra course at Polytechnique for this one ❤️
Vector operator+(const Vector &a, const Vector &b)
{
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b)
{
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(double scalar, const Vector &vector)
{
    return Vector(scalar * vector[0], scalar * vector[1], scalar * vector[2]);
}
Vector operator*(const Vector &vector, double scalar)
{
    return Vector(vector[0] * scalar, vector[1] * scalar, vector[2] * scalar);
}
Vector operator*(const Vector &a, const Vector &b)
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector &vector, double scalar)
{
    return Vector(vector[0] / scalar, vector[1] / scalar, vector[2] / scalar);
}

double dot(const Vector &a, const Vector &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector &a, const Vector &b)
{
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}