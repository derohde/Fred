cmake_minimum_required(VERSION 3.16)
project(Fred LANGUAGES CXX C)

set( CMAKE_CXX_STANDARD 14 )

if(NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
if(NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive") #supress error in older gcc
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-trapping-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-loop")

include_directories(${CMAKE_SOURCE_DIR}/include)

add_definitions(-D_GLIBCXX_PARALLEL)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
    add_compile_definitions(WITH_OMP)
else()
    message(WARNING "Compiling without openmp")
endif()

add_subdirectory(pybind11)

pybind11_add_module(backend
            src/curve.cpp
            src/point.cpp
            src/interval.cpp
            src/frechet.cpp
            src/jl_transform.cpp
            src/simplification.cpp
            src/dynamic_time_warping.cpp
            src/clustering.cpp
            src/config.cpp
            src/bounding.cpp
            src/fred_python_wrapper.cpp
)
