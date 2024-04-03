cmake_minimum_required(VERSION 3.10)
project(Raytracer)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add source files
set(SOURCES
    main.cpp
)

# Add header files
include_directories(stb)

# Add executable
add_executable(Raytracer ${SOURCES})
