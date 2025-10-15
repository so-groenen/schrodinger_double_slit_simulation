cmake_minimum_required(VERSION 3.10.0)
project(schrodinger VERSION 0.1.0 LANGUAGES C CXX)

add_executable(schrodinger src/schrodinger_equation.cpp
                           src/main.cpp)
find_package(eigen3 REQUIRED)
find_package(raylib REQUIRED)

if(eigen3_FOUND)
    message(STATUS "Found eigen3")
else()
    message(STATUS "Did NOT Found eigen3")
endif()
if(raylib_FOUND)
    message(STATUS "Found raylib")
else()
    message(STATUS "Did NOT Found raylib")
endif()

target_include_directories(${PROJECT_NAME} PRIVATE include)
target_compile_features   (${PROJECT_NAME} PRIVATE cxx_std_23)
target_compile_options    (${PROJECT_NAME} PRIVATE -O3 -march=native)
target_link_libraries     (${PROJECT_NAME} PRIVATE raylib stdc++exp)