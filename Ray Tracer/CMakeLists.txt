cmake_minimum_required(VERSION 3.10)
project(Raytracer)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add source files
set(SOURCES
    main.cpp
)

# Add header files
include_directories(stb)

# Add folder cat_files containing cat.obj
file(COPY cat_files DESTINATION ${CMAKE_BINARY_DIR})

# Add executable
add_executable(Raytracer ${SOURCES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Raytracer PUBLIC OpenMP::OpenMP_CXX)
endif()
