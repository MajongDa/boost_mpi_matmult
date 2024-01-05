#include "mpi_matmult.hpp"

using namespace std;

// вывод матрицы на экран
void matrix_out ( double* A, int N1,int N2 )
    {
    for ( int i = 0; i < N1; i++ )
        {
        for ( int j = 0; j < N2; j++ )
            cout << A[i*N2+j] << ' ';
        cout << endl;
        }
    }


int main ( int argc, char* argv [] )
{
    int n = atoi(argv [ 1 ]);

    float *Ap, *Bp, *Cp, s;
    vec Av ( n * n ), Bv ( n * n ), Cv ( n * n );
    matr Am ( n, n ), Bm ( n, n ), Cm ( n, n );
    Ap = ( float * ) calloc ( n * n, sizeof ( float) );
    Bp = ( float* ) calloc ( n * n, sizeof ( float ) );
    Cp = ( float * ) calloc ( n * n, sizeof ( float ) );
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < n; j++ )
            {
            Ap [i * n + j] = Av(i * n + j) = Am(i, j) =
            Bp [i * n + j] = Bv(i * n + j) = Bm(i, j) = ( i == j ? 4 : 1 );
            Cp [i * n + j] = 0;
            }

    mpi::environment env ( argc, argv );

    matmult ( Ap, Bp, Cp, n );
    free(Ap); free(Bp); free(Cp);

    matmult( Av, Bv, Cv, n);

    matmult(Am, Bm, Cm, n);
}
