#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "nrutil.h"

/*
Skopiowane z ksiazki Numerical Recipes
*/
//void choldc(float **a, int n, float p[])
//{
//    void nrerror(char error_text[]);
//    int i,j,k;
//    float sum;
//
//    for (i=1;i<=n;i++)
//    {
//        for (j=i;j<=n;j++)
//        {
//            for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
//            if (i == j)
//            {
//                if (sum <= 0.0) nrerror("choldc failed");
//                p[i]=sqrt(sum);
//            } else a[j][i]=sum/p[i];
//        }
//    }
//}
//void cholsl(float **a, int n, float p[], float b[], float x[])
//{
//    int i,k;
//    float sum;
//
//    for (i=1;i<=n;i++)
//    {
//        for (sum=b[i],k=i-1;k>=1;k--) sum -= a[i][k]*x[k];
//        x[i]=sum/p[i];
//    }
//    for (i=n;i>=1;i--)
//    {
//        for (sum=x[i],k=i+1;k<=n;k++) sum -= a[k][i]*x[k];
//        x[i]=sum/p[i];
//    }
//}

/*
Wersja pierwsza ze slajd�w
*/
void choldc(double**a, double**l, int n)
{
    void nrerror(char error_text[]);
    int i,j,k;
    double sum;
    clock_t begin, end;
    double time_spent;

    begin = clock();
    for (k = 1; k <= n; k++)
    {
        for (sum = a[k][k], j = 1; j <= k - 1; j++) sum -= l[k][j] * l[k][j];
        l[k][k] = sqrt(sum);
        for (i = k + 1; i <= n; i++)
        {
            for (sum = a[i][k], j = 1; j <= k - 1; j++) sum -= l[i][j] * l[k][j];
            if (k == i)
            {
                if (sum <= 0.0) nrerror("choldc failed");
            } else l[i][k] = sum / l[k][k];
        }
    }
    end = clock();
    printf("time: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
}

/*
Wersja druga ze slajd�w
*/
void choldc2(double**a, double**l, int n)
{
    void nrerror(char error_text[]);
    int i,j,k;
    double sum;
    clock_t begin, end;
    double time_spent;

    begin = clock();
    for (k = 1; k <= n - 1; k++)
    {
        l[k][k] = sqrt(a[k][k]);
        for (i = k + 1; i <= n; i++) l[i][k] = a[i][k] / l[k][k];
        for (j = k + 1; j <= n; j++)
        {
            for (i = j; i <= n; i++) a[i][j] = a[i][j] - l[i][k] * l[j][k];
        }
    }
    l[n][n] = sqrt(a[n][n]);
    end = clock();
    printf("time: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
}




int main(void)
{
//    int n = 3;
//    double**l;
//    l = matrix(1, 3, 1, 3);
//    double**a;
//    a = matrix(1, 3, 1, 3);
//    //przyklad z angielskiej wiki http://en.wikipedia.org/wiki/Cholesky_decomposition#Example
//    a[1][1] = 4;
//    a[2][1] = 12;
//    a[3][1] = -16;
//    a[1][2] = 12;
//    a[2][2] = 37;
//    a[3][2] = -43;
//    a[1][3] = -16;
//    a[2][3] = -43;
//    a[3][3] = 98;
//
//    choldc(a, l, n);
//    printf("%f\t0\t\t0\n", l[1][1], l[1][2], l[1][3]);
//    printf("%f\t%f\t0\n", l[2][1], l[2][2], l[2][3]);
//    printf("%f\t%f\t%f\n", l[3][1], l[3][2], l[3][3]);
//
//    choldc2(a, l, n);
//    printf("\n\n%f\t0\t\t0\n", l[1][1], l[1][2], l[1][3]);
//    printf("%f\t%f\t0\n", l[2][1], l[2][2], l[2][3]);
//    printf("%f\t%f\t%f\n", l[3][1], l[3][2], l[3][3]);

    //przyklad macierzy 1000x1000 wygenerowanej w matlabie
    FILE * pFile;
    pFile = fopen ("SPDmatrix1000","r");
    FILE * pFile2;
    pFile2 = fopen ("Lower1000","r");
    FILE * pFile3;
    pFile3 = fopen ("result","w");

    double f;
    int n = 1000;
    double**a;
    double**l1, **l2, **l_reference;
    a = dmatrix(1, 1000, 1, 1000);
    l1 = dmatrix(1, 1000, 1, 1000);
    l2 = dmatrix(1, 1000, 1, 1000);
    l_reference = dmatrix(1, 1000, 1, 1000);
    int i, j;
    //wczytujemy macierz A do faktoryzacji
    //i macierz L wyliczana w matlabie
    for (i = 1; i <= 1000; i++)
    {
        for (j = 1; j <= 1000; j++)
        {
            fscanf(pFile, "%lf", &f);
            a[i][j] = f;
            fscanf(pFile2, "%lf", &f);
            l_reference[i][j] = f;
        }
    }

    choldc(a, l1, n);
    //liczymy blad jako srednie odchylenie
    //dla elementow pod diagonala
    double error = 0;
    int counter = 0;
    for (i = 1; i <= 1000; i++)
    {
        for (j = 1; j <= 1000; j++)
        {
            if (i < j) continue;
            counter++;
            error += (l1[i][j] - l_reference[i][j]) * (l1[i][j] - l_reference[i][j]);
            fprintf(pFile3, "%lf\t", l1[i][j]);
        }
        fprintf(pFile3, "\n");
    }
    printf("error: %lf\n", sqrt(error / (double)counter));

    choldc2(a, l2, n);
    error = 0;
    counter = 0;
    for (i = 1; i <= 1000; i++)
    {
        for (j = 1; j <= 1000; j++)
        {
            if (i < j) continue;
            counter++;
            error += (l2[i][j] - l_reference[i][j]) * (l2[i][j] - l_reference[i][j]);
        }
    }
    printf("error: %lf\n", sqrt(error / (double)counter));

return 0;
}
