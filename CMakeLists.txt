cmake_minimum_required(VERSION 3.21)
project(mpi_matmult_test LANGUAGES CXX)

find_package(MPI REQUIRED)
find_package(Boost 1.82.0 REQUIRED)


include_directories(${Boost_INCLUDE_DIR})

add_compile_options(-O3)

add_library ( mpi_matmult mpi_matmult.cpp )
target_link_libraries ( mpi_matmult boost_mpi MPI::MPI_CXX)

add_executable ( mpi_matmult_test mpi_matmult_test.cpp )
target_link_libraries( mpi_matmult_test mpi_matmult boost_mpi MPI::MPI_CXX)


