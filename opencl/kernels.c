#define A(x,y) A[x*dimension + y]
#define L(x,y) L[x*dimension + y]


#pragma OPENCL EXTENSION cl_amd_printf : enable

__kernel void choldc_gpu(__global float* A, __global float* L, const unsigned int dimension)
{
int x = get_global_id(0);
int i, k, j;
float sum;

    //zeby workery ktore wykraczaja poza obszar macierzy nic nie kombinowaly
    if (x < dimension)
    {
//        //to robia wszyscy sekwencyjnie
        for (k = 0; k < dimension; k++)
        {
            //tu sie nic nie zrownolegli
            sum = A(k,k);
            for (j = 0; j < k - 1; j++) sum = sum - L(k,j) * L(k,j);
            L(k,k) = sqrt(sum);
            //za to tutaj mozna uzyc workerow wiekszych od aktualnego k
            //do przeliczenia rownolegle wszystkich wierszy pod elementem na diagonali
            if (x >= k)
            {
                //kazdy worker bierze po jednym i ( po jednym wierszu )
                i = x + 1;
                sum = A(i,k);
                for (j = 0; j < k - 1; j++) sum = sum - L(i,j) * L(k,j);
                L(i,k) = sum / L(k,k);
            }
        }
    }

}
