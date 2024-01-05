#include "mpi_matmult.hpp"

template void matmult<int> ( int *, int *, int *, int);
template void matmult<float> ( float *, float*, float *, int);
template void matmult<double> ( double*, double *, double *, int);
template <class T>
void matmult( T* A, T* B, T* C, int n )
    {


    mpi::communicator world;

    int p = world.size();
    int r = world.rank();

    int cnt = n / p;
    int from = r * cnt;
    int to = n;


    mpi::timer t;
    t.restart();

    if ( r != p-1 ) to = from + cnt;

    for ( int i = from; i < to; ++i )
        for ( int j = 0; j < n;  ++j )
            for ( int k = 0; k < n; ++k )
               C[ i * n + j ] += A[ i * n + k ] * B[ k * n + j ];

    world.barrier ();

    if ( r != 0 )
        {
        world.send ( 0, 1, from );
        world.send ( 0, 2, to );
        world.send ( 0, 3, C + from * n, n * ( to - from ) );
        //MPI_Send ( &in, 1, MPI_INT, 0, 1, MPI_COMM_WORLD );
        //MPI_Send ( &ik,1,MPI_INT,0,2,MPI_COMM_WORLD );
        //MPI_Send ( C+in*n, ( ik-in ) *n,MPI_DOUBLE,0,0,MPI_COMM_WORLD );
        }
    else
        {
        for ( int i = 1; i < p; ++i )
            {
            world.recv ( i, 1, from );
            world.recv ( i, 2, to );
            world.recv ( i, 3, C + from * n, n * ( to - from ) );
            //MPI_Recv ( &in, 1, MPI_INT, ii, 1, MPI_COMM_WORLD, &st );
            //MPI_Recv ( &ik,1,MPI_INT,ii,2,MPI_COMM_WORLD, &st );
            //MPI_Recv ( C+in*n, ( ik-in ) *n, MPI_DOUBLE,ii,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE );
            }

        std::cout<<"Итоговое время счёта с * = " << t.elapsed() << std::endl;
        /*matrix_out ( A, n, n );
        cout << "b= " << endl;
        matrix_out ( B, n, n );
        cout<<"c= "<<endl;
        matrix_out ( C, n, n );*/
        }

}


void matmult( vec A, vec B, vec C, int n )
    {

    mpi::communicator world;

    int p = world.size();
    int r = world.rank();

    int cnt = n / p;
    int from = r * cnt;
    int to = n;


    mpi::timer t;
    t.restart();

    if ( r != p-1 ) to = from + cnt;

    for ( int i = from; i < to; ++i )
        for ( int j = 0; j < n;  ++j )
            for ( int k = 0; k < n; ++k )
               C[ i * n + j ] += A[ i * n + k ] * B[ k * n + j ];

    world.barrier ();

    if ( r != 0 )
        {
        world.send ( 0, 1, from );
        world.send ( 0, 2, to );
        world.send ( 0, 3, &C( from * n ), ( to - from) * n );
        }
    else
        {
        for ( int i = 1; i < p; ++i )
            {
            world.recv ( i, 1, from );
            world.recv ( i, 2, to );
            world.recv ( i, 3, &C( from * n ), ( to - from) * n );
            }

        std::cout<<"Итоговое время счёта с uBLAS::vector = " << t.elapsed() << std::endl;
        /*for ( auto i: C)
            cout <<  i << "\t" ;
        cout << endl;*/
        }
    }


void matmult( matr A, matr B, matr C, int n )
    {

    mpi::communicator world;

    int p = world.size();
    int r = world.rank();

    int cnt = n / p;
    int from = r * cnt;
    int to = n;


    mpi::timer t;
    t.restart();

    if ( r != p-1 ) to = from + cnt;

    for ( int i = from; i < to; ++i )
        for ( int j = 0; j < n;  ++j )
            for ( int k = 0; k < n; ++k )
               C( i, j ) += A( i, k ) * B( k, j );

    world.barrier ();

    if ( r != 0 )
        {
        world.send ( 0, 1, from );
        world.send ( 0, 2, to );
        world.send ( 0, 3, &C( from ), ( to - from ) * n );
        }
    else
        {
        for ( int i = 1; i < p; ++i )
            {
            world.recv ( i, 1, from );
            world.recv ( i, 2, to );
            world.recv ( i, 3, &C( from ), ( to - from ) * n );
            }

        std::cout<<"Итоговое время счёта с uBLAS::matrix = " << t.elapsed() << std::endl;
        /*for ( auto i: C)
            cout <<  i << "\t" ;
        cout << endl;*/
        }
    }
