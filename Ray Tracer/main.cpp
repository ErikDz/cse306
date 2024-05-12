#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <chrono>
#include <iostream>
#include <cfloat>
#include <tuple>
#include <list>
//#include "omp.h"

/*
NOTE TO SELF:
REMEMBER TO CD TO BUILD BEFORE BUILDING!!!
*/

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






class Geometry {
public:
    virtual Intersection intersect(const Ray& ray) = 0;
};


class BoundingBox {
public:
    BoundingBox() {};
    Vector minPoint;
    Vector maxPoint;

    explicit BoundingBox(Vector min, Vector max) {
        minPoint = min;
        maxPoint = max;
    }
};

struct Node {
    BoundingBox bounding_box;
    Node* c_left;
    Node* c_right;

    int start_trig;
    int end_trig;

    Node() : c_left(nullptr), c_right(nullptr), start_trig(0), end_trig(0) {}
};

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

class TriangleMesh : public Geometry {
    double scaling_factor;
    Vector translation;
    Vector color;
    double refractive_index;
    bool reflects;

public:

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    Node* root;
    BoundingBox bounding_box;

    ~TriangleMesh() {}
    TriangleMesh(double scaling_factor,
                 Vector translation,
                 Vector color = Vector(0., 0., 0.),
                 double refractive_index = 1.,
                 bool reflects = false) {
        this->scaling_factor = scaling_factor;
        this->translation = translation;
        this->color = color;
        this->refractive_index = refractive_index;
        this->reflects = reflects;
        this->root = new Node;
    };


    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);

                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }

            }

        }
        fclose(f);
        this-> BuildBVH(root, 0, indices.size());
    }



// AUX FOR THE BOUNDING BOX
BoundingBox calculateBoundingBox() {
    double minX{DBL_MAX}, minY{DBL_MAX}, minZ{DBL_MAX};
    double maxX{DBL_MIN}, maxY{DBL_MIN}, maxZ{DBL_MIN};
    for (const auto& vertex : vertices) {
        Vector scaledVertex = scaling_factor * vertex + translation;
        minX = std::min(minX, scaledVertex[0]);
        maxX = std::max(maxX, scaledVertex[0]);
        minY = std::min(minY, scaledVertex[1]);
        maxY = std::max(maxY, scaledVertex[1]);
        minZ = std::min(minZ, scaledVertex[2]);
        maxZ = std::max(maxZ, scaledVertex[2]);
    }

    return BoundingBox(Vector(minX, minY, minZ), Vector(maxX, maxY, maxZ));
}

std::tuple<double, double> calculateRayPlaneIntersection(const Vector& normal, const Vector& min, const Vector& max, const Ray& ray) {
    double tMin = dot(min - ray.origin, normal) / dot(ray.direction, normal);
    double tMax = dot(max - ray.origin, normal) / dot(ray.direction, normal);
    double t0 = std::min(tMin, tMax);
    double t1 = std::max(tMin, tMax);
    return std::make_tuple(t0, t1);
}

bool doesRayIntersectBoundingBox(const Ray& ray, BoundingBox bounding_box, double* t) {
    BoundingBox box = calculateBoundingBox();

    // Check for ray-plane intersection along each axis
    auto [txMin, txMax] = calculateRayPlaneIntersection(Vector(1, 0, 0), box.minPoint, box.maxPoint, ray);
    auto [tyMin, tyMax] = calculateRayPlaneIntersection(Vector(0, 1, 0), box.minPoint, box.maxPoint, ray);
    auto [tzMin, tzMax] = calculateRayPlaneIntersection(Vector(0, 0, 1), box.minPoint, box.maxPoint, ray);

    double tMin = std::min(std::min(txMax, tyMax), tzMax);
    double tMax = std::max(std::max(txMin, tyMin), tzMin);

    *t = (tMin> tMax) ? tMax : *t;
    return (tMin > tMax);
}

Vector CalculateBarycenter(int triangleIndex) {
    // Calculate the vertices of the triangle with the given index
    Vector v1 = (scaling_factor * this->vertices[this->indices[triangleIndex].vtxi]) + translation;
    Vector v2 = (scaling_factor * this->vertices[this->indices[triangleIndex].vtxj]) + translation;
    Vector v3 = (scaling_factor * this->vertices[this->indices[triangleIndex].vtxk]) + translation;
    
    // Return the barycenter (average) of the three vertices
    return (v1 + v2 + v3) / 3.0;
}

void PartitionTriangles(int startTriangle, int endTriangle, int& pivotIndex, int longestAxis, const Vector& middlePoint) {
    // Partition the triangles based on the barycenter along the longest axis
    for (int i = startTriangle; i < endTriangle; i++) {
        Vector barycenter = this->CalculateBarycenter(i);
        if (barycenter[longestAxis] < middlePoint[longestAxis]) {
            std::swap(indices[i], indices[pivotIndex]);
            pivotIndex++;
        }
    }
}

void BuildBVH(Node* node, int startTriangle, int endTriangle) {
    // Set the bounding box and triangle range for the current node
    node->bounding_box = this->calculateBoundingBox();
    node->start_trig = startTriangle;
    node->end_trig = endTriangle;

    // Calculate the diagonal and middle point of the bounding box
    Vector diagonal = node->bounding_box.maxPoint - node->bounding_box.minPoint;
    Vector middlePoint = node->bounding_box.minPoint + diagonal * 0.5;
    
    // Determine the longest axis of the bounding box
    int longestAxis = (abs(diagonal[0]) > abs(diagonal[1]) &&
                       abs(diagonal[0]) > abs(diagonal[2])) ? 0 :
                      (abs(diagonal[1]) > abs(diagonal[2])) ? 1 : 2;

    int pivotIndex = startTriangle;

    // Partition the triangles
    PartitionTriangles(startTriangle, endTriangle, pivotIndex, longestAxis, middlePoint);

    // Check if the partition is valid or if the triangle range is too small
    if (pivotIndex <= startTriangle
        || pivotIndex >= endTriangle - 1
        || endTriangle - startTriangle < 5)
        return;

    // Create child nodes and recursively build the BVH for each child
    node->c_left = new Node;
    node->c_right = new Node;
    this->BuildBVH(node->c_left, startTriangle, pivotIndex);
    this->BuildBVH(node->c_right, pivotIndex, endTriangle);
}


// INTERSECTION CODE ================================================
Intersection intersect(const Ray& ray) override {
    //if (!doesRayIntersectBoundingBox(ray))
        //return Intersection();

    Intersection intersection(this->color, this->refractive_index, this->reflects);
    //Vector vertexA, vertexB, vertexC, edge1, edge2, normal;
    double t;
    double tMin{DBL_MAX};
    if (!doesRayIntersectBoundingBox(ray, this->bounding_box, &t))
        return Intersection();

    std::list<Node*> nodesToVisit;
    nodesToVisit.push_front(root);

    // we start the while loop
    while (!nodesToVisit.empty()){
        Node *currentNode = nodesToVisit.back();
        nodesToVisit.pop_back();

        // if the current node is a leaf 
        if (currentNode->c_left){
            // if the ray intersects the bounding box of the current node
            if (doesRayIntersectBoundingBox(ray, currentNode->c_left->bounding_box, &t)){
                if (t<tMin){
                    nodesToVisit.push_back(currentNode->c_left);
                }
            }
            // same for the right child
            if (doesRayIntersectBoundingBox(ray, currentNode->c_right->bounding_box, &t)){
                if (t<tMin){
                    nodesToVisit.push_back(currentNode->c_right);
                }
            }
        } else{
            // else

            Vector vertexA, vertexB, vertexC, edge1, edge2, normal;
            for (const auto& index : indices) {
                vertexA = scaling_factor * vertices[index.vtxi] + translation;
                vertexB = scaling_factor * vertices[index.vtxj] + translation;
                vertexC = scaling_factor * vertices[index.vtxk] + translation;
                edge1 = vertexB - vertexA;
                edge2 = vertexC - vertexA;
                normal = cross(edge1, edge2);

                double beta = dot(cross(vertexA - ray.origin, ray.direction), edge2) / dot(ray.direction, normal);
                double gamma = -dot(cross(vertexA - ray.origin, ray.direction), edge1) / dot(ray.direction, normal);
                double alpha = 1.0 - beta - gamma;
                double t = dot(vertexA - ray.origin, normal) / dot(ray.direction, normal);

                if (alpha >= 0 && beta >= 0 && gamma >= 0 && t > 0 && t < tMin) {
                    tMin = t;
                    intersection.intersects = true;
                    intersection.distance = t;
                    intersection.position = vertexA + beta * edge1 + gamma * edge2;
                    intersection.normal = normal;
                }
            }
        }
    }
    return intersection;
}
};


// Sphere class representing a sphere with its properties
class Sphere : public Geometry {
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
    Intersection intersect(const Ray& ray) override{
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
    std::vector<Geometry*> geometries;
    Vector lightSource;
    double lightIntensity = 1e5;
    //double lightIntensity = 1e4;

public:
    explicit Scene(Vector lightSource) { this->lightSource = lightSource; }
    void addGeometry(Geometry* geometry) { geometries.push_back(geometry); }

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
        for (auto& sphere : geometries) {
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
            //Ray randomRay = Ray(localPosition, randomCosineDirection(localNormal));
            //color = color + intersection.color * getColor(randomRay, depth - 1);
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

    //timer
    auto start = std::chrono::high_resolution_clock::now();

    // Create a scene with a light source
    Scene scene = Scene(Vector(-10, 20, 40));

    // Create spheres with different properties
//    Sphere* mirror = new Sphere(Vector(20, 0, 0), 10, Vector(1., 1., 1.), true, 1.5);
//    Sphere* refracted = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
    Sphere* hollowOuter = new Sphere(Vector(0, 0, 0), 10, Vector(1., 1., 1.), false, 1.5);
    Sphere* hollowInner = new Sphere(Vector(0, 0, 0), 9.8, Vector(1., 1., 1.), false, 1.5, true);
//    Sphere* solid = new Sphere(Vector(-20, 0, 0), 9.8, Vector(1., 1., 1.), false, 0, false);
//
//    scene.addGeometry(mirror);
//    scene.addGeometry(refracted);
    scene.addGeometry(hollowOuter);
    scene.addGeometry(hollowInner);
//

    // Create spheres for the room
    Sphere* ceiling = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
    Sphere* floor = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
    Sphere* front = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere* back = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
    Sphere* left = new Sphere(Vector(1000, 0, 0), 940, Vector(0, 1, 1));
    Sphere* right = new Sphere(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));

    scene.addGeometry(ceiling);
    scene.addGeometry(floor);
    scene.addGeometry(front);
    scene.addGeometry(back);
    scene.addGeometry(left);
    scene.addGeometry(right);
    
    //Sphere* solid = new Sphere(Vector(0, 0, 0), 9.8, Vector(1, 1, 1));
    //scene.addGeometry(solid);
    //TriangleMesh* cat = new TriangleMesh(0.6, Vector(0, -10, 0), Vector(1., 1., 1.), 1.0, false);
    //cat.readOBJ("cat_files/cat.obj");
    //cat.readOBJ("cat_files/cat.obj");
    //cat->readOBJ("cat_files/cat.obj");
    //scene.addGeometry(cat);

    //int imageWidth = 1024;
    //int imageHeight = 1024;
    int imageWidth = 1024;
    int imageHeight = 1024;
    std::vector<unsigned char> image(imageWidth * imageHeight * 3, 0);
    Vector camera = Vector(0, 0, 55);
    double fieldOfView = degToRad(60);
    double gamma = 2.2;
    int maxDepth = 20;
    int raysPerPixel = 600;

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
            //pixelColor = pixelColor/raysPerPixel;

            image[(y * imageWidth + x) * 3 + 0] = std::min(255., pow(pixelColor[0] / raysPerPixel, 1. / gamma) * 255);
            image[(y * imageWidth + x) * 3 + 1] = std::min(255., pow(pixelColor[1] / raysPerPixel, 1. / gamma) * 255);
            image[(y * imageWidth + x) * 3 + 2] = std::min(255., pow(pixelColor[2] / raysPerPixel, 1. / gamma) * 255);
        }

    // stop timing!
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time taken by function: " << duration.count() << " milliseconds" << std::endl;


    stbi_write_png("image.png", imageWidth, imageHeight, 3, &image[0], 0);
    return 0;
}