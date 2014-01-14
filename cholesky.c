#include "cholesky.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>

/*
Wersja pierwsza ze slajdów
*/
double** choldc(double **A, double **L, int dimension)
{
    int i,j,k;
    double sum;
    clock_t begin, end;

omp_set_dynamic(0);
omp_set_num_threads(4);

    begin = clock();
    for (k = 1; k <= dimension; k++)
    {
        sum = A[k][k];
        for (j = 1; j <= k - 1; j++) sum = sum - (L[k][j] * L[k][j]);
        L[k][k] = sqrt(sum);
        #pragma omp parallel for private(i, j, sum) shared(A, L, dimension)
        for (i = k + 1; i <= dimension; i++)
        {
            sum = A[i][k];
            for (j = 1; j <= k - 1; j++) sum = sum - L[i][j] * L[k][j];
            L[i][k] = sum / L[k][k];
        }
    }
    end = clock();
    printf("Time method 1: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);

return L;
}

/*
Wersja druga ze slajdów
*/
double** choldc2(double **A, double **L, int dimension)
{
    int i,j,k;
    clock_t begin, end;

omp_set_dynamic(0);
omp_set_num_threads(16);

    begin = clock();
    for (k = 1; k <= dimension - 1; k++)
    {
        L[k][k] = sqrt(A[k][k]);
        #pragma omp parallel for private(i) shared(L, A, dimension)
        for (i = k + 1; i <= dimension; i++) L[i][k] = A[i][k] / L[k][k];
        #pragma omp parallel for private(j, i) shared(L, A, dimension)
        for (j = k + 1; j <= dimension; j++)
        {
            for (i = j; i <= dimension; i++) A[i][j] = A[i][j] - L[i][k] * L[j][k];
        }
    }
    L[dimension][dimension] = sqrt(A[dimension][dimension]);
    end = clock();
    printf("Time method 2: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);

return L;
}


