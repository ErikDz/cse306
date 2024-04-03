#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <chrono>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "./stb/stb_image.h"
#define _CRT_SECURE_NO_WARNINGS 1

// Random engine and distribution for generating random numbers
static std::default_random_engine randomEngine(10);
static std::uniform_real_distribution<double> randomDistribution(0.0, 1.0);

// Function to convert degrees to radians
double degToRad(double degrees) {
    return degrees * M_PI / 180.0;
}

// Vector class representing a 3D vector
class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        components[0] = x;
        components[1] = y;
        components[2] = z;
    }
    double squaredNorm() const {
        return components[0] * components[0] + components[1] * components[1] + components[2] * components[2];
    }
    double norm() const {
        return sqrt(squaredNorm());
    }
    Vector normalized() {
        double n = norm();
        Vector normalizedVector;
        normalizedVector[0] = components[0] / n;
        normalizedVector[1] = components[1] / n;
        normalizedVector[2] = components[2] / n;
        return normalizedVector;
    }
    double operator[](int i) const { return components[i]; };
    double& operator[](int i) { return components[i]; };
    double components[3];
};

// Vector arithmetic operations
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(double scalar, const Vector& vector) {
    return Vector(scalar * vector[0], scalar * vector[1], scalar * vector[2]);
}
Vector operator*(const Vector& vector, double scalar) {
    return Vector(vector[0] * scalar, vector[1] * scalar, vector[2] * scalar);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector& vector, double scalar) {
    return Vector(vector[0] / scalar, vector[1] / scalar, vector[2] / scalar);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

// Intersection structure representing the intersection point and its properties
struct Intersection {
    Vector position;
    Vector normal;
    Vector color;
    double distance = 0.;
    double refractiveIndex;
    bool isReflective;
    bool intersects = false;

    Intersection(Vector color = Vector(0., 0., 0.),
                 double refractiveIndex = 1.,
                 bool isReflective = false) {
        this->color = color;
        this->refractiveIndex = refractiveIndex;
        this->isReflective = isReflective;
    }
};

// Ray class representing a ray with an origin and direction
class Ray {
public:
    Vector origin;
    Vector direction;
    explicit Ray(Vector origin, Vector direction) {
        this->origin = origin;
        this->direction = direction;
    }
};

// Sphere class representing a sphere with its properties
class Sphere {
private:
    Vector center;
    Vector color;
    double radius;
    double refractiveIndex;
    bool isReflective;
    bool isHollow;

public:
    explicit Sphere(Vector center,
                    double radius,
                    Vector color,
                    bool isReflective = false,
                    double refractiveIndex = 1.,
                    bool isHollow = false) {
        this->center = center;
        this->radius = radius;
        this->color = color;
        this->isReflective = isReflective;
        this->refractiveIndex = refractiveIndex;
        this->isHollow = isHollow;
    }

    // Method to compute the intersection of a ray with the sphere
    Intersection intersect(const Ray& ray) {
        Intersection intersection(this->color, this->refractiveIndex, this->isReflective);

        Vector originToCenter = ray.origin - this->center;
        double discriminant = pow(dot(ray.direction, originToCenter), 2)
                              - (dot(originToCenter, originToCenter) - pow(radius, 2));

        if (discriminant >= 0.) {
            double t1 = dot(ray.direction, -1. * originToCenter) - sqrt(discriminant);
            double t2 = dot(ray.direction, -1. * originToCenter) + sqrt(discriminant);
            intersection.distance = (t1 > 0) ? t1 : ((t2 > 0) ? t2 : 0.0);
            intersection.intersects = (t2 < 0.) ? false : true;
        }

        intersection.position = ray.origin + (intersection.distance * ray.direction);
        intersection.normal = (intersection.position - this->center).normalized();
        intersection.normal = (this->isHollow) ? -1. * intersection.normal : intersection.normal;
        return intersection;
    }
};

// Scene class representing the scene with spheres and light source
class Scene {
private:
    std::vector<Sphere*> spheres;
    Vector lightSource;
    double lightIntensity = 1e5;

public:
    explicit Scene(Vector lightSource) { this->lightSource = lightSource; }
    void addSphere(Sphere* sphere) { spheres.push_back(sphere); }

    // Method to generate a random direction based on the cosine-weighted distribution
    Vector randomCosineDirection(const Vector& normal) {
        double r1 = randomDistribution(randomEngine);
        double r2 = randomDistribution(randomEngine);
        double x = sqrt(1 - r2) * cos(2. * M_PI * r1);
        double y = sqrt(1 - r2) * sin(2. * M_PI * r1);
        double z = sqrt(r2);

        double minValue = std::numeric_limits<double>::max();
        int axis = 0;
        for (int i = 0; i < 3; i++)
            if (abs(normal[i]) < minValue) {
                minValue = abs(normal[i]);
                axis = i;
            }

        Vector tangent1 = (axis == 0) ? Vector(0., normal[2], -normal[1]).normalized()
                                      : (axis == 1) ? Vector(normal[2], 0., -normal[0]).normalized()
                                                    : Vector(normal[1], -normal[0], 0.).normalized();
        Vector tangent2 = cross(normal, tangent1);
        return tangent1 * x + tangent2 * y + normal * z;
    }

    // Method to compute the intersection of a ray with the scene
    Intersection intersect(const Ray& ray) {
        Intersection intersectionMain, intersectionTemp;
        double minDistance = std::numeric_limits<double>::max();
        for (auto& sphere : spheres) {
            intersectionTemp = sphere->intersect(ray);
            if (intersectionTemp.intersects && intersectionTemp.distance < minDistance) {
                minDistance = intersectionTemp.distance;
                intersectionMain = intersectionTemp;
            }
        }
        return intersectionMain;
    }


    // Method to compute the color of a ray by recursive ray tracing
    Vector getColor(const Ray& ray, int depth) {
        if (depth < 0) return Vector(0., 0., 0.);

        Intersection intersection = intersect(ray);
        Vector color(0., 0., 0.);

        if (intersection.intersects) {
            double epsilon = 1e-10;
            Vector localPosition = intersection.position + (epsilon * intersection.normal);
            Vector localNormal = intersection.normal;

            if (intersection.isReflective) {
                Ray reflectedRay = Ray(localPosition, ray.direction - (2 * dot(ray.direction, localNormal) * localNormal));
                return getColor(reflectedRay, depth - 1);
            }

            if (intersection.refractiveIndex != 1.) {
                double Nu = dot(ray.direction, localNormal);
                double n1 = (Nu > 0.) ? intersection.refractiveIndex : 1.;
                double n2 = (Nu > 0.) ? 1. : intersection.refractiveIndex;
                localNormal = (Nu > 0.) ? -1. * localNormal : localNormal;

                localPosition = intersection.position - (epsilon * localNormal);
                Nu = dot(ray.direction, localNormal);
                if (1. - pow(n1 / n2, 2) * (1. - pow(Nu, 2)) > 0.) {
                    // Fresnel's law
                    Vector transmissionDirection = (n1 / n2) * (ray.direction - Nu * localNormal);
                    Vector normalDirection = -1. * localNormal * sqrt(1. - pow(n1 / n2, 2) * (1 - pow(Nu, 2)));
                    Vector refractedDirection = transmissionDirection + normalDirection;
                    double k0 = pow((n1 - n2) / (n1 + n2), 2);
                    double reflectionProbability = k0 + (1 - k0) * pow(1 - abs(dot(localNormal, refractedDirection)), 5);
                    if (randomDistribution(randomEngine) < reflectionProbability) {
                        Ray reflectedRay = Ray(localPosition, ray.direction - (2 * dot(ray.direction, intersection.normal) * intersection.normal));
                        return getColor(reflectedRay, depth - 1);
                    }
                    else {
                        Ray refractedRay = Ray(localPosition, refractedDirection);
                        return getColor(refractedRay, depth - 1);
                    }
                }
                else {
                    //total internal reflection
                    Ray reflectedRay = Ray(localPosition, ray.direction - (2 * dot(intersection.normal, ray.direction) * intersection.normal));
                    return getColor(reflectedRay, depth - 1);
                }
            }

            // Add direct lighting in diffuse case
            double distance = (lightSource - localPosition).norm();
            Vector lightDirection = (lightSource - localPosition).normalized();
            Intersection lightIntersection = intersect(Ray(lightSource, lightDirection * (-1.)));
            double visibility = (!lightIntersection.intersects || lightIntersection.distance > distance) ? 1. : 0.;
            color = lightIntensity / (4 * M_PI * distance * distance) * intersection.color / M_PI * visibility * std::max(0., dot(lightDirection, localNormal));

            // Indirect lighting
            Ray randomRay = Ray(localPosition, randomCosineDirection(localNormal));
            color = color + intersection.color * getColor(randomRay, depth - 1);
        }

        return color;
    }
};

// Function to generate random numbers using the Box-Muller transform
void BoxMuller(double stdev, double& x, double& y) {
    double r1 = randomDistribution(randomEngine);
    double r2 = randomDistribution(randomEngine);
    x = stdev * sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
    y = stdev * sqrt(-2 * log(r1)) * sin(2 * M_PI * r2);
}

int main() {
    // Create a scene with a light source
    Scene scene = Scene(Vector(-10, 20, 40));

    // Create spheres with different properties
    Sphere* mirror = new Sphere(Vector(20, 0, 0), 10, Vector(1., 1., 1.), true, 1.5);
    Sphere* refracted = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
    Sphere* hollowOuter = new Sphere(Vector(-20, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
    Sphere* hollowInner = new Sphere(Vector(-20, 0, 0), 9.8, Vector(1., 1., 1.), false, 1.5, true);
    // Create spheres for the room
    Sphere* ceiling = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
    Sphere* floor = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
    Sphere* front = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere* back = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
    Sphere* left = new Sphere(Vector(1000, 0, 0), 940, Vector(0, 1, 1));
    Sphere* right = new Sphere(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));

    // Add spheres to the scene
    scene.addSphere(mirror);
    scene.addSphere(refracted);
    scene.addSphere(hollowOuter);
    scene.addSphere(hollowInner);
    scene.addSphere(ceiling);
    scene.addSphere(floor);
    scene.addSphere(front);
    scene.addSphere(back);
    scene.addSphere(left);
    scene.addSphere(right);

    int imageWidth = 512;
    int imageHeight = 512;
    std::vector<unsigned char> image(imageWidth * imageHeight * 3, 0);
    Vector camera = Vector(0, 0, 55);
    double fieldOfView = degToRad(60);
    double gamma = 2.2;
    int maxDepth = 5;
    int raysPerPixel = 10;

    #pragma omp parallel for schedule(dynamic, 1)
    for (int y = 0; y < imageHeight; y++)
        for (int x = 0; x < imageWidth; x++) {
            Vector pixelColor = Vector(0., 0., 0.);
            double offsetX, offsetY;

            for (int i = 0; i < raysPerPixel; i++) {
                BoxMuller(0.5, offsetX, offsetY);
                Vector pixel = Vector(camera[0] + (x + offsetX) + 0.5 - imageWidth / 2,
                                      camera[1] - (y + offsetY) - 0.5 + imageHeight / 2,
                                      camera[2] - imageWidth / (2 * tan(fieldOfView / 2)));
                Ray ray = Ray(camera, (pixel - camera).normalized());
                pixelColor = pixelColor + scene.getColor(ray, maxDepth);
            }

            image[(y * imageWidth + x) * 3 + 0] = std::min(255., pow(pixelColor[0] / raysPerPixel, 1. / gamma) * 255);
            image[(y * imageWidth + x) * 3 + 1] = std::min(255., pow(pixelColor[1] / raysPerPixel, 1. / gamma) * 255);
            image[(y * imageWidth + x) * 3 + 2] = std::min(255., pow(pixelColor[2] / raysPerPixel, 1. / gamma) * 255);
        }

    stbi_write_png("image.png", imageWidth, imageHeight, 3, &image[0], 0);
    return 0;
}