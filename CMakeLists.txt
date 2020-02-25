cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(Fred LANGUAGES CXX C)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -shared")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive") #supress error in older gcc
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-trapping-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-loop")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-undefined")

include_directories(${CMAKE_SOURCE_DIR}/include)

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Boost 1.63 COMPONENTS system chrono ${BPY} ${BNPY} REQUIRED)
find_package(OpenMP REQUIRED)

add_definitions(-D_GLIBCXX_PARALLEL)

include_directories(${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

if(OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
endif()
if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

endif()
link_libraries(OpenMP::OpenMP_CXX)

PYTHON_ADD_MODULE(backend
            src/fred_python_wrapper.cpp
            src/curve.cpp
            src/point.cpp
            src/interval.cpp
            src/frechet.cpp
)
