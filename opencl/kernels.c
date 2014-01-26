#define A(x,y) A[x*dimension + y]
#define L(x,y) L[x*dimension + y]

__kernel void choldc_gpu(__global float* A, __global float* L, const unsigned int dimension)
{
int i = get_global_id(0);
int k, j;
float sum;

if (i < dimension)
{
    for (k = 0; k < dimension; k++)
    {
        //to robimy po kolei - tylko pierwszy watek
        if (i == 0)
        {
            sum = A(k,k);
            for (j = 0; j < k - 1; j++) sum = sum - L(k,j) * L(k,j);
            L(k,k) = sqrt(sum);
        }
    }
    //wszystkie workery czekaja na zerowego zeby wyliczyl pierwiastek
    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    for (k = 0; k < dimension; k++)
    {
        if (i > k)
        {
            sum = A(i,k);
            for (j = 0; j < k - 1; j++) sum = sum - L(i,j) * L(k,j);
            L(i,k) = sum / L(k,k);
        }
    }
}
}
