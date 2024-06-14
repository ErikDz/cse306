cmake_minimum_required(VERSION 3.10)

# Project name
project(td9)

set(SOURCES
    td9.cpp
    classes.cpp
)

set(HEADERS
    classes.h
)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add executable
add_executable(td9 td9.cpp)

# Additional compiler flags (optional)
if (CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(td9 PRIVATE -Wall -Wextra -pedantic)
endif()
