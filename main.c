#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <emmintrin.h>
#include "nrutil.h"
#include "spd_matrix.h"
#include "cholesky.h"


int main(int argc, char* argv[])
{
    srand( time( NULL ) );
	int i,j;
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
    L = choldc_openMP(A, L, dimension);
    L_t = clone_matrix(L, dimension);
    L_t = transpose_matrix(L_t, dimension);

    //Oblicza norme nowej macierzy
    A_clone = multiply(L, L_t, A_clone, dimension);
    norm2 = frobenius_norm(A_clone, dimension);

    //Faktoryzacja Choleskyego metoda 2.
    L = choldc2_openMP(A, L, dimension);
    L_t = clone_matrix(L, dimension);
    L_t = transpose_matrix(L_t, dimension);

    //Oblicza norme nowej macierzy
    A_clone = multiply(L, L_t, A_clone, dimension);
    norm3 = frobenius_norm(A_clone, dimension);


    printf("Norm1: % 20.16lf\n", norm1);
    printf("Norm2: % 20.16lf\n", norm2);
    printf("Norm3: % 20.16lf\n", norm3);
    printf("Error method 1: % 20.16lf\n", fabs(norm1 - norm2));
    printf("Error method 2: % 20.16lf\n", fabs(norm1 - norm3));


//    //alokujemy pamiec z wyrownaniem
//    __m128d *data, *L;
//    data = _mm_malloc(dimension * dimension * sizeof(double), 16);
//    L = _mm_malloc(dimension * dimension * sizeof(double), 16);
//    //ladujemy tam macierz w zwektoryzowanej postaci
//    int count = convert_to_sse(data, A, dimension);
//
////    double *temp = (double*) data;
////    for (i = 0; i < count * 2; i++)
////        printf("% 20.16lf\n", temp[i]);
//
//    choldc_sse(data, L, A, dimension);

return 0;
}
