#ifndef SPD_MATRIX
#define SPD_MATRIX

double **dmatrix();
double random_double(double fMin, double fMax);
void print_matrix(double** A, int dimension);
double** generate_random_matrix(double** A, int dimension);
double** clone_matrix(double** A, int dimension);
double** transpose_matrix(double** A, int dimension);
double** construct_symetric_matrix(double** A, int dimension);
double** create_identity_matrix(int dimension);
double** matrix_positive_definite(double** A, int dimension);
double** create_lower_triangular(double **A, int dimension);
double** multiply(double **L, double **L_t, double **A, int dimension);
double frobenius_norm(double** L, int dimension);

#endif // SPD_MATRIX
