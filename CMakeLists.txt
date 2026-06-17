cmake_minimum_required(VERSION 3.10)
project(OBC_MPI)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(MPI REQUIRED)
find_package(Eigen3 REQUIRED) # Décommente si besoin

# Cible pointant vers le dossier src/
add_executable(OBC_MPI src/main.cpp)
target_link_libraries(OBC_MPI PRIVATE Eigen3::Eigen MPI::MPI_CXX)
