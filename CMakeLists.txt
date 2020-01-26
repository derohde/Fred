cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(Fred LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -shared")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive") #supress error in older gcc
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-undefined")

include_directories(${CMAKE_SOURCE_DIR}/include)

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Boost 1.63 COMPONENTS chrono ${BPY} ${BNPY} REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES} ${Boost_LIBRARIES} OpenMP::OpenMP_CXX)

add_definitions(-D_GLIBCXX_PARALLEL)

PYTHON_ADD_MODULE(Fred
            src/fred_python_wrapper.cpp
            src/curve.cpp
            src/point.cpp
            src/interval.cpp
            src/frechet.cpp
)
