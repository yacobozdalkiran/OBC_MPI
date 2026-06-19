cmake_minimum_required(VERSION 3.10)
project(OBC_MPI)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()


find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED) 
include(FetchContent)
FetchContent_Declare(
        lime
        GIT_REPOSITORY https://github.com/usqcd-software/c-lime.git
        GIT_TAG        master  
)
FetchContent_MakeAvailable(lime)

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
    src/flow/gradient_flow.cpp
    src/io/params.h
    src/io/io.cpp
    src/io/ildg.cpp
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

add_executable(test_ecmc src/bin/test_ecmc.cpp
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
    src/io/params.h
    src/ecmc/ecmc.cpp
    src/flow/gradient_flow.cpp
    src/io/io.cpp
    src/io/ildg.cpp
)
target_link_libraries(test_ecmc PRIVATE MPI::MPI_CXX Eigen3::Eigen lime OpenMP::OpenMP_CXX)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(test_ecmc PRIVATE
            -O3
            -xHost            
            -qopenmp
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(test_ecmc PRIVATE
            -O3
            #-ffast-math
            -march=native
            -Wall -Wextra -Wpedantic
            -fopenmp
    )
endif()
 
add_executable(test_hb src/bin/test_hb.cpp
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
    src/io/params.h
    src/heatbath/heatbath.cpp
    src/heatbath/heatbath_mpi.cpp
    src/flow/gradient_flow.cpp
    src/io/io.cpp
    src/io/ildg.cpp
)
target_link_libraries(test_hb PRIVATE MPI::MPI_CXX Eigen3::Eigen lime OpenMP::OpenMP_CXX)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(test_hb PRIVATE
            -O3
            -xHost            
            -qopenmp
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(test_hb PRIVATE
            -O3
            #-ffast-math
            -march=native
            -Wall -Wextra -Wpedantic
            -fopenmp
    )
endif()

add_executable(test_ildg src/bin/test_ildg.cpp
    src/su3/utils.cpp
    src/gauge/GaugeField.cpp
    src/geometry/Geometry.cpp
    src/geometry/types.h
    src/mpi/MpiTopology.cpp
    src/mpi/HalosExchange.cpp
    src/mpi/HalosShift.h
    src/mpi/Shift.cpp
    src/io/ildg.cpp
)
target_link_libraries(test_ildg PRIVATE MPI::MPI_CXX Eigen3::Eigen lime OpenMP::OpenMP_CXX)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(test_ildg PRIVATE
            -O3
            -xHost            
            -qopenmp
            -fp-model precise
            -Wall -Wextra
    )
else()
    target_compile_options(test_ildg PRIVATE
            -O3
            -march=native
            -Wall -Wextra -Wpedantic
            -fopenmp
    )
endif()
