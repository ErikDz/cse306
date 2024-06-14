cmake_minimum_required(VERSION 3.10)
project(Geometry)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add source files
set(SOURCES
    main.cpp
    utils/classes.cpp
    voronoi.cpp
    fluid.cpp
    rendering.cpp
    gallouet.cpp
    utils/svg.cpp
)

# Add header files
include_directories(stb)
include_directories(liblbfgs)
include_directories(utils)
include_directories(Eigen)

# Add the header files for proper compilation
set(HEADERS
    voronoi.h
    fluid.h
    rendering.h
    gallouet.h
    utils/classes.h
    utils/svg.h
)

# Add executable
add_executable(Geometry ${SOURCES} ${HEADERS})

# Link the lbfgs library if it has a separate implementation file
add_library(lbfgs liblbfgs/lbfgs.c)
target_link_libraries(Geometry lbfgs)