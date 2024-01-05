#ifndef MPI_MATMULT_HPP_INCLUDED
#define MPI_MATMULT_HPP_INCLUDED

#include <boost/mpi.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>


namespace mpi = boost::mpi;
using matr = boost::numeric::ublas::matrix<float>;
using vec = boost::numeric::ublas::vector<float>;

template<typename T>
void matmult( T* , T* , T* , int );

void matmult(vec, vec, vec, int );

void matmult( matr , matr , matr , int );

#endif // MPI_MATMULT_HPP_INCLUDED
