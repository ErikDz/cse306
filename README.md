# README - Ray Tracing Renderer

This project is a ray tracing renderer implemented in C++. It utilizes the principles of ray tracing to generate realistic 3D images by simulating the interaction of light with various objects in a scene.

## Dependencies

- C++ compiler with C++11 support
- OpenMP library for parallel rendering
- `stb_image_write.h` and `stb_image.h` header files from the stb library for image I/O

## Usage

1. Clone the repository and navigate to the project directory.
2. Compile the code using a C++ compiler with C++11 support and OpenMP enabled.
3. Run the compiled executable.
4. The rendered image will be saved as `image.png` in the project directory.

## Acknowledgments

- The stb library (https://github.com/nothings/stb) for providing the `stb_image_write.h` and `stb_image.h` header files used for image I/O.
- The OpenMP library (https://www.openmp.org/) for enabling parallel rendering.

---

cd build
cmake ..
cmake --build .