#define _CRT_SECURE_NO_WARNINGS 1
#include <algorithm>
#include <random>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "./stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "./stb_image.h"

constexpr int kNumIterations = 100;


struct Vector {
    double x, y, z;

    Vector(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    double SquaredNorm() const { return x * x + y * y + z * z; }
    double Norm() const { return sqrt(SquaredNorm()); }
    Vector Normalize() const {
        double norm = Norm();
        return Vector(x / norm, y / norm, z / norm);
    }
};

Vector operator+(const Vector& a, const Vector& b) { return Vector(a.x + b.x, a.y + b.y, a.z + b.z); }
Vector operator-(const Vector& a, const Vector& b) { return Vector(a.x - b.x, a.y - b.y, a.z - b.z); }
Vector operator*(double scalar, const Vector& v) { return Vector(scalar * v.x, scalar * v.y, scalar * v.z); }
Vector operator*(const Vector& v, double scalar) { return scalar * v; }
Vector operator*(const Vector& a, const Vector& b) { return Vector(a.x * b.x, a.y * b.y, a.z * b.z); }
Vector operator/(const Vector& v, double scalar) { return Vector(v.x / scalar, v.y / scalar, v.z / scalar); }
double Dot(const Vector& a, const Vector& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vector Cross(const Vector& a, const Vector& b) {
    return Vector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}


Vector RandomDirection() {
    double r1 = static_cast<double>(rand()) / RAND_MAX;
    double r2 = static_cast<double>(rand()) / RAND_MAX;
    double x = cos(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double y = sin(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double z = 1 - 2 * r2;
    return Vector(x, y, z);
}

void ColorMatching(int width, int height, int num_channels, int num_color_channels,
                   unsigned char* image_data, unsigned char* color_data) {
    size_t num_pixels = width * height;
    std::vector<std::pair<double, int>> image_projections(num_pixels);
    std::vector<std::pair<double, int>> color_projections(num_pixels);

    for (size_t i = 0; i < kNumIterations; i++) {
        Vector random_vector = RandomDirection();

        for (size_t j = 0; j < num_pixels; j++) {
            Vector image_pixel(image_data[j * num_channels], image_data[j * num_channels + 1], image_data[j * num_channels + 2]);
            Vector color_pixel(color_data[j * num_color_channels], color_data[j * num_color_channels + 1], color_data[j * num_color_channels + 2]);
            image_projections[j] = std::make_pair(Dot(image_pixel, random_vector), j);
            color_projections[j] = std::make_pair(Dot(color_pixel, random_vector), j);
        }

        std::sort(image_projections.begin(), image_projections.end());
        std::sort(color_projections.begin(), color_projections.end());

        for (size_t j = 0; j < num_pixels; j++) {
            int index = image_projections[j].second;
            Vector image_pixel(image_data[index * num_channels], image_data[index * num_channels + 1], image_data[index * num_channels + 2]);
            Vector advected_pixel = image_pixel + (color_projections[j].first - image_projections[j].first) * random_vector;
            image_data[index * num_channels] = static_cast<unsigned char>(advected_pixel.x);
            image_data[index * num_channels + 1] = static_cast<unsigned char>(advected_pixel.y);
            image_data[index * num_channels + 2] = static_cast<unsigned char>(advected_pixel.z);
        }
    }
}

/*  ------------------------ MAIN ---------------------------------   */

int main(int argc, char** argv) {
    int image_width, image_height, image_channels;
    int color_width, color_height, color_channels;

    unsigned char* image_data = stbi_load("imgA.jpg", &image_width, &image_height, &image_channels, 0);
    if (!image_data) {
        // Handle image loading error
        return 1;
    }

    unsigned char* color_data = stbi_load("redim.jpg", &color_width, &color_height, &color_channels, 0);
    if (!color_data) {
        // Handle color source loading error
        stbi_image_free(image_data);
        return 1;
    }

    ColorMatching(image_width, image_height, image_channels, color_channels, image_data, color_data);

    stbi_write_png("output.png", image_width, image_height, image_channels, image_data, 0);

    stbi_image_free(image_data);
    stbi_image_free(color_data);

    return 0;
}
