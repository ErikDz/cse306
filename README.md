# README - Ray Tracing Renderer

This project is a ray tracing renderer implemented in C++. It utilizes the principles of ray tracing to generate realistic 3D images by simulating the interaction of light with various objects in a scene.

## Features

- Supports rendering of spheres with different properties such as reflectivity, refraction, and hollowness.
- Implements direct lighting and indirect lighting for accurate shading.
- Supports reflection and refraction of light rays.
- Uses the Box-Muller transform for generating random numbers.
- Utilizes the cosine-weighted distribution for generating random directions.
- Implements Fresnel's law for determining the reflection probability.
- Supports rendering of a room with walls, ceiling, and floor.
- Utilizes OpenMP for parallel rendering to improve performance.

## Dependencies

- C++ compiler with C++11 support
- OpenMP library for parallel rendering
- `stb_image_write.h` and `stb_image.h` header files from the stb library for image I/O

## Usage

1. Clone the repository and navigate to the project directory.
2. Compile the code using a C++ compiler with C++11 support and OpenMP enabled.
3. Run the compiled executable.
4. The rendered image will be saved as `image.png` in the project directory.

## Possible References

Here are some references that one could have used while writing this ray tracing renderer:

1. "Physically Based Rendering: From Theory to Implementation" by Matt Pharr, Wenzel Jakob, and Greg Humphreys
   - This book provides a comprehensive guide to physically based rendering techniques, including ray tracing, shading, and global illumination.

2. "Ray Tracing in One Weekend" by Peter Shirley
   - This book offers a beginner-friendly introduction to ray tracing and guides the reader through the implementation of a basic ray tracer in a weekend.

3. "Realistic Ray Tracing" by Peter Shirley and R. Keith Morley
   - This book delves into advanced ray tracing techniques, including Monte Carlo methods, importance sampling, and advanced shading models.

4. "Fundamentals of Computer Graphics" by Peter Shirley and Steve Marschner
   - This book covers the fundamental concepts of computer graphics, including coordinate systems, transformations, and shading models.

5. "Ray Tracing Gems" edited by Eric Haines and Tomas Akenine-Möller
   - This book is a collection of articles and techniques related to ray tracing, covering topics such as acceleration structures, sampling, and optimization.

6. "Raytracing.github.io" (https://raytracing.github.io/)
   - This website provides a series of articles and tutorials on ray tracing, covering topics from basic ray tracing to advanced techniques.

7. "Scratchapixel" (https://www.scratchapixel.com/)
   - This website offers tutorials and articles on various computer graphics topics, including ray tracing, shading, and global illumination.

These references can provide valuable insights and techniques for implementing a ray tracing renderer and understanding the underlying concepts.

## Acknowledgments

- The stb library (https://github.com/nothings/stb) for providing the `stb_image_write.h` and `stb_image.h` header files used for image I/O.
- The OpenMP library (https://www.openmp.org/) for enabling parallel rendering.

Feel free to explore and modify the code to experiment with different scenes, materials, and rendering techniques. Happy ray tracing!
