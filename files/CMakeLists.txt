cmake_minimum_required(VERSION 3.14)

# Project name and version
project(HybridDocking VERSION 1.0
	    DESCRIPTION "Docking with MPI and OpenMP"
		LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Add the executables
add_executable(serial_docking src/serial_docking.cpp)
add_executable(hybrid_docking src/hybrid_docking.cpp)

# Add include for serial_docking
target_include_directories(serial_docking PRIVATE include/)
target_include_directories(hybrid_docking PRIVATE include/)

# Find OpenMP and MPI
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
target_link_libraries(hybrid_docking PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX)
