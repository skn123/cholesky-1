#include "cholesky.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <emmintrin.h>
#include "nrutil.h"

/*
Wersja pierwsza ze slajdów
*/
double** choldc(double **A, double **L, int dimension)
{
    void nrerror(char error_text[]);
    int i,j,k;
    double sum;
    clock_t begin, end;
    double time_spent;

    begin = clock();
    for (k = 1; k <= dimension; k++)
    {
        for (sum = A[k][k], j = 1; j <= k - 1; j++) sum -= L[k][j] * L[k][j];
        L[k][k] = sqrt(sum);
        for (i = k + 1; i <= dimension; i++)
        {
            for (sum = A[i][k], j = 1; j <= k - 1; j++) sum -= L[i][j] * L[k][j];
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
    double time_spent;

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
    printf("Time method 2: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);

return L;
}

void choldc_sse(__m128d *data, __m128d *L, int dimension)
{
    int i,j,k, count = 0;
    clock_t begin, end;
    double time_spent;
    __m128d temp1, temp2;
    __m128d *L_sse = (__m128d*) L;
    __m128d *data_sse = data;


    begin = clock();
    for (k = 1; k <= dimension - 1; k += 2)
    {
        temp1 = _mm_sqrt_pd(*data_sse);
        memcpy(L_sse, &temp1, 128);
        data_sse++;
        L_sse++;
        for (i = k + 1; i <= dimension; i++)
        {
            temp2 = _mm_div_pd(*data_sse, temp1);
            memcpy(L_sse, &temp2, 128);
            data_sse++;
            L_sse++;
        }
    }
    end = clock();
    printf("Time method SSE: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);

    double *kij = (double*) L;
    for (i = 0; i < 12; i += 2)
        printf("% 20.16lf\t\t% 20.16lf\n", kij[i], kij[i + 1]);
}
