#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

void choldc(float **a, int n, float p[])
{
    void nrerror(char error_text[]);
    int i,j,k;
    float sum;

    for (i=1;i<=n;i++)
    {
        for (j=i;j<=n;j++)
        {
            for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
            if (i == j)
            {
                if (sum <= 0.0) nrerror("choldc failed");
                p[i]=sqrt(sum);
            } else a[j][i]=sum/p[i];
        }
    }
}

void cholsl(float **a, int n, float p[], float b[], float x[])
{
    int i,k;
    float sum;

    for (i=1;i<=n;i++)
    {
        for (sum=b[i],k=i-1;k>=1;k--) sum -= a[i][k]*x[k];
        x[i]=sum/p[i];
    }
    for (i=n;i>=1;i--)
    {
        for (sum=x[i],k=i+1;k<=n;k++) sum -= a[k][i]*x[k];
        x[i]=sum/p[i];
    }
}


int main(void)
{
    int n = 3;
    float *x;
    x = vector(1, 3);
    float *p;
    p = vector(1, 3);
    float *b;
    b = vector(1, 3);
    b[1] = 2;
    b[2] = 8;
    b[3] = 10;
    float **a;
    a = matrix(1, 3, 1, 3);
    a[1][1] = 1;
    a[2][1] = 2;
    a[3][1] = 3;
    a[1][2] = 2;
    a[2][2] = 8;
    a[3][2] = 10;
    a[1][3] = 3;
    a[2][3] = 10;
    a[3][3] = 22;

    choldc(a, n, p);
    cholsl(a, n, p, b, x);

    printf("%f %f %f\n", x[1], x[2], x[3]);

    return 0;
}
