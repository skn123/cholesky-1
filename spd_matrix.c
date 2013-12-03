#include "spd_matrix.h"
#include "nrutil.h"
#include <math.h>
#include <stdlib.h>
#include <emmintrin.h>

double random_double(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void print_matrix(double** A, int dimension)
{
	int i,j;

	for(i = 1; i <= dimension; i++)
    {
    	for(j = 1; j <= dimension; j++)
    	{
    		printf("%20.16lf\t\t\t",A[i][j]);
    	}
      	printf("\n");
    }
}

double** generate_random_matrix(double** A, int dimension)
{
	int i,j;

	for (i = 1; i <= dimension; i++)
    {
        for (j = 1; j <= dimension; j++)
        {
            A[i][j] = random_double(0, 1);
            //printf("Result: %lf\n", A[i][j]);
        }
    }
	return A;
}

double** clone_matrix(double** A, int dimension)
{
	double** cloned_matrix = dmatrix(1, dimension, 1, dimension);
	int i,j;
	for(i = 1; i <= dimension; i++)
    {
    	for(j = 1; j <= dimension; j++)
    	{
        	cloned_matrix[i][j] = A[i][j];
        	//printf("%f, ",cloned_matrix[i][j]);
        }
    }
    return cloned_matrix;
}

double** transpose_matrix(double** A, int dimension)
{
	int i,j;
	double** A_t;
	A_t = clone_matrix(A, dimension);
	for(i = 1; i <= dimension; i++) //transponowanie macierzy
    {
    	for(j = 1; j <= dimension; j++)
        	A[j][i] = A_t[i][j];
    }
    return A;
}

double** construct_symetric_matrix(double** A, int dimension)
{
	double** temp;
	temp = clone_matrix(A, dimension);
	double** A_t;
	int i,j;
	A_t = transpose_matrix(temp, dimension);
	//print_matrix(A_t, dimension);
	for(i = 1; i <= dimension; i++) // A = A+A'
    	for(j = 1; j <= dimension; j++)
      		A[i][j] = A[i][j] + A_t[i][j];

    return A;
}

double** create_identity_matrix(int dimension) //n*I(n)
{
	double** nI = dmatrix(1, dimension, 1, dimension);
	int i,j;
	for(i = 1; i <= dimension; i++)
	{
    	for(j = 1; j <= dimension; j++)
    	{
    		if (i == j)
      			nI[i][j] = dimension;
      		else
      			nI[i][j] = 0;
      	}
    }
    return nI;
}

double** matrix_positive_definite(double** A, int dimension)
{
	double** nI;
	int i,j;
	nI = create_identity_matrix(dimension);

	//print_matrix(I,dimension);
	for(i = 1; i <= dimension; i++) // A = A + n*I(n);
    	for(j = 1; j <= dimension; j++)
      		A[i][j] = A[i][j] + nI[i][j];

    return A;
}

double** create_lower_triangular(double **A, int dimension)
{
    double p, temp;
    int i, j, k;

    for(i = 1; i <= dimension; i++)
    {
        for(j = 2; j <= dimension - i; j++)
        {
            if(A[i + j][i] / A[i][i])
            {
                p = A[i + j][i] / A[i][i];
                for(k = 1; k <= dimension - i; k++)
                {
                    temp = A[i][i + k] * p;
                    A[i + j][i + k] -= temp;
                }
            }
        }
    }
    return A;
}

double** multiply(double **L, double **L_t, double **A, int dimension)
{
    int i, j, k;

    for (i = 1; i <= dimension; i++)
    {
        for(j = 1; j <= dimension; j++)
        {
            A[i][j] = 0;
            for(k = 1; k <= dimension; k++)
                A[i][j] += L[i][k] * L_t[k][j];
        }
    }

return A;
}

double frobenius_norm(double** L, int dimension)
{
    int i, j;
    double norm = 0.0;

    for (i = 1; i <= dimension; i++)
    {
        for (j = 1; j <= dimension; j++)
            norm += L[i][j] * L[i][j];
    }

return sqrt(norm);
}


int convert_to_sse(__m128d *data, double **A, int dimension)
{
    int i, j, count = 0;
    double *temp;
    temp = malloc(128);

//    data = _mm_malloc(128, 16);
//    __m128d kij1;
//    data = memcpy(data, &temp[1], 128);
//    kij1 = _mm_add_pd(*data, *data);
//    temp = memcpy(temp, &kij1, 128);
//
//    printf("a[0]: % 20.16lf\na[1]: % 20.16lf", temp[0], temp[1]);


    //ladujemy tam elementy macierzy A w odpowiedniej kolejnosci
    //lykamy po dwie kolumny
    for (i = 1; i <= dimension; i += 2)
    {
        temp[0] = A[i][i];
        temp[1] = A[i + 1][i + 1];
        //kopiujemy dwa double
        memcpy(data + count, temp, 128);
        count++;
        for (j = i + 1; j <= dimension; j++)
        {
            //bierzemy po dwie kolumny ale pomiewaz
            //mamy macierz trojkatna to druga kolumna ma zawsze
            //jeden element mniej, stad 0 jako padding
            if (j == dimension)
            {
                temp[0] = A[j][i];
                temp[1] = 0.0;
                memcpy(data + count, temp, 128);
                count++;
            }
            else
            {
                temp[0] = A[j][i];
                temp[1] = A[j + 1][i + 1];
                memcpy(data + count, temp, 128);
                count++;
            }
        }
    }

return count;
}











