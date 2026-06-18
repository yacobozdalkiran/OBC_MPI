cmake_minimum_required(VERSION 3.10)
project(OBC_MPI)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED) 

add_executable(test_mpi src/bin/test_mpi.cpp
    src/su3/utils.cpp
    src/gauge/GaugeField.cpp
    src/geometry/Geometry.cpp
    src/geometry/types.h
    src/observables/observables.cpp
    src/observables/observables_mpi.cpp
    src/mpi/MpiTopology.cpp
    src/mpi/HalosExchange.cpp
    src/mpi/HalosShift.h
    src/mpi/Shift.cpp
)
target_link_libraries(test_mpi PRIVATE MPI::MPI_CXX Eigen3::Eigen lime OpenMP::OpenMP_CXX)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(test_mpi PRIVATE
            -O3
            -xHost            
            -qopenmp
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(test_mpi PRIVATE
            -O3
            #-ffast-math
            -march=native
            -Wall -Wextra -Wpedantic
            -fopenmp
    )
endif()
