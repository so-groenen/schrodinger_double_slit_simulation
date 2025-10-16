cmake_minimum_required(VERSION 3.10.0)
project(double_slit_experiment VERSION 0.1.0 LANGUAGES C CXX)

add_executable(${PROJECT_NAME}  src/crank_nicolson_builder.cpp
                                src/gaussian_wavefunction_builder.cpp
                                src/schrodinger_equation_builder.cpp
                                src/schrodinger_equation.cpp
                                src/interferometer.cpp
                                src/main.cpp)

find_package(eigen3)
if(eigen3_FOUND)
    message(STATUS "Found eigen3")
else()
    message(STATUS "Did NOT Found eigen3: Fetchting from gitlab")
    include(FetchContent)
    FetchContent_Declare(
    eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG d426838d011d7f938ab0c077dd48b7bf1bb71f59
    )
    FetchContent_MakeAvailable(eigen3)
endif()

find_package(raylib)
if(raylib_FOUND)
    message(STATUS "Found raylib")
else()
    message(STATUS "Did NOT Found raylib: Fetchting from github")
    include(FetchContent)
    FetchContent_Declare(
    raylib
    GIT_REPOSITORY https://github.com/raysan5/raylib.git
    GIT_TAG c1ab645ca298a2801097931d1079b10ff7eb9df8 
    )
    FetchContent_MakeAvailable(eigen3)
endif()

target_include_directories(${PROJECT_NAME} PRIVATE include)
target_compile_features   (${PROJECT_NAME} PRIVATE cxx_std_23)
target_link_libraries     (${PROJECT_NAME} PRIVATE raylib Eigen3::Eigen)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options    (${PROJECT_NAME} PRIVATE /O2 /arch:AVX2)
else()
    target_compile_options    (${PROJECT_NAME} PRIVATE -O3 -march=native)
    target_link_libraries     (${PROJECT_NAME} PRIVATE stdc++exp)
endif()