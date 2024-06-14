cmake_minimum_required(VERSION 3.10)

# Project name
project(td9)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add source files
set(SOURCES
    td9.cpp
    classes.cpp
)

# Add header files (optional, not required for compilation)
set(HEADERS
    classes.h
)

# Add executable
add_executable(td9 ${SOURCES})

# Additional compiler flags (optional)
if (CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(td9 PRIVATE -Wall -Wextra -pedantic)
endif()
