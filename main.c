#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include "spd_matrix.h"
#include "cholesky.h"


int main(int argc, char* argv[])
{
    srand( time( NULL ) );
	if (!argv[1])
	{
	    printf("Specify matrix dimension.\n");
	    exit(-1);
	}
	int dimension = atoi(argv[1]);
    printf("Dimension: %d\n", dimension);
    double norm1, norm2, norm3;
    double** A = dmatrix(1, dimension, 1, dimension);
    double** A_clone = dmatrix(1, dimension, 1, dimension);
    double** L = dmatrix(1, dimension, 1, dimension);
    double** L_t = dmatrix(1, dimension, 1, dimension);

    //Generowanie macierzy SPD
	A = generate_random_matrix(A, dimension);
    A = construct_symetric_matrix(A, dimension);
	A = matrix_positive_definite(A, dimension);
    norm1 = frobenius_norm(A, dimension);
//    print_matrix(A, dimension);

    //Faktoryzacja Choleskyego metoda 1.
    L = choldc(A, L, dimension);
    L_t = clone_matrix(L, dimension);
    L_t = transpose_matrix(L_t, dimension);

    //Oblicza norme nowej macierzy
    A_clone = multiply(L, L_t, A_clone, dimension);
    norm2 = frobenius_norm(A_clone, dimension);

    //Faktoryzacja Choleskyego metoda 2.
    L = choldc2(A, L, dimension);
    L_t = clone_matrix(L, dimension);
    L_t = transpose_matrix(L_t, dimension);

    //Oblicza norme nowej macierzy
    A_clone = multiply(L, L_t, A_clone, dimension);
    norm3 = frobenius_norm(A_clone, dimension);

    printf("Error method 1: % 20.16lf\n", fabs(norm1 - norm2));
    printf("Error method 2: % 20.16lf\n", fabs(norm1 - norm3));


return 0;
}
