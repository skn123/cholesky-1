#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

/*
Skopiowane z ksiazki Numerica Recepies
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
Wersja pierwsza ze slajdów
*/
void choldc(float **a, float **l, int n)
{
    void nrerror(char error_text[]);
    int i,j,k;
    float sum;

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
}

/*
Wersja druga ze slajdów
*/
void choldc2(float **a, float **l, int n)
{
    void nrerror(char error_text[]);
    int i,j,k;
    float sum;

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
}




int main(void)
{
    int n = 3;
    float **l;
    l = matrix(1, 3, 1, 3);
    float **a;
    a = matrix(1, 3, 1, 3);
    //przyklad z angielskiej wiki http://en.wikipedia.org/wiki/Cholesky_decomposition#Example
    a[1][1] = 4;
    a[2][1] = 12;
    a[3][1] = -16;
    a[1][2] = 12;
    a[2][2] = 37;
    a[3][2] = -43;
    a[1][3] = -16;
    a[2][3] = -43;
    a[3][3] = 98;

    choldc(a, l, n);
    printf("%f\t0\t\t0\n", l[1][1], l[1][2], l[1][3]);
    printf("%f\t%f\t0\n", l[2][1], l[2][2], l[2][3]);
    printf("%f\t%f\t%f\n", l[3][1], l[3][2], l[3][3]);

    choldc2(a, l, n);
    printf("\n\n%f\t0\t\t0\n", l[1][1], l[1][2], l[1][3]);
    printf("%f\t%f\t0\n", l[2][1], l[2][2], l[2][3]);
    printf("%f\t%f\t%f\n", l[3][1], l[3][2], l[3][3]);

return 0;
}
