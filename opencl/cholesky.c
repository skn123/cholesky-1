#include "cholesky.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/*
Wersja pierwsza ze slajdów
*/
float** choldc(float **A, float **L, int dimension)
{
    int i,j,k;
    float sum;
    clock_t begin, end;

    begin = clock();
    for (k = 1; k <= dimension; k++)
    {
        sum = A[k][k];
        for (j = 1; j <= k - 1; j++) sum = sum - (L[k][j] * L[k][j]);
        L[k][k] = sqrt(sum);
        for (i = k + 1; i <= dimension; i++)
        {
            sum = A[i][k];
            for (j = 1; j <= k - 1; j++) sum = sum - L[i][j] * L[k][j];
            L[i][k] = sum / L[k][k];
        }
    }
    end = clock();
    printf("Time \tCPU: % 20.16lf\n", (float)(end - begin) / CLOCKS_PER_SEC);

return L;
}

/*
Wersja druga ze slajdów
*/
float** choldc2(float **A, float **L, int dimension)
{
    int i,j,k;
    clock_t begin, end;

    begin = clock();
    for (k = 1; k <= dimension - 1; k++)
    {
        L[k][k] = sqrt(A[k][k]);
        for (i = k + 1; i <= dimension; i++) L[i][k] = A[i][k] / L[k][k];
        for (j = k + 1; j <= dimension; j++)
        {
            for (i = j; i <= dimension; i++) A[i][j] = A[i][j] - L[i][k] * L[j][k];
        }
    }
    L[dimension][dimension] = sqrt(A[dimension][dimension]);
    end = clock();
    printf("Time \tCPU: % 20.16lf\n", (float)(end - begin) / CLOCKS_PER_SEC);

return L;
}


