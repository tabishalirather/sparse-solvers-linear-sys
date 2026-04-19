#include <stdio.h>
#include <stdlib.h>
#include "../Library/csr.h"

void matrix_mult_csr(){
printf("Callimg csr_matrix_mult\n");
}

int main() {
    printf("printing in sheet1\n");
    //part a: Allocating memory for matrix.

    enum dims { ROWS = 3, COLS = 5 };

int* mat_A = malloc(ROWS * COLS * sizeof(int));

// row-major order
// 0  3  1  0  0   → indices 0–4
// 1  0  4  0  5   → indices 5–9
// 0  3  2  4  0   → indices 10–14

printf("malloc for matA returns: %p\n", mat_A);

// Row 0
mat_A[0 * COLS + 0] = 0;
mat_A[0 * COLS + 1] = 3;
mat_A[0 * COLS + 2] = 1;
mat_A[0 * COLS + 3] = 0;
mat_A[0 * COLS + 4] = 0;

// Row 1
mat_A[1 * COLS + 0] = 1;
mat_A[1 * COLS + 1] = 0;
mat_A[1 * COLS + 2] = 4;
mat_A[1 * COLS + 3] = 0;
mat_A[1 * COLS + 4] = 5;

// Row 2
mat_A[2 * COLS + 0] = 0;
mat_A[2 * COLS + 1] = 3;
mat_A[2 * COLS + 2] = 2;
mat_A[2 * COLS + 3] = 4;
mat_A[2 * COLS + 4] = 0;
for (int i=0; i<ROWS; i++) {
    for (int j=0; j<COLS; j++) {
        printf("%d ", *(mat_A + i*COLS +j));
    }
    printf("\n");
}
free(mat_A);
typedef struct csr_matrix
{
  int nr;          /*!< Number of rows. */
  int nc;          /*!< Number of columns. */
  int nz;          /*!< Number of non-zero entries. */
  int *row;        /*!< Pointers to the beginning of each row. */
  int *col;        /*!< Column indices of non-zero coefficients. */
  double *coeff;   /*!< Non-zero coefficients. */
} csr;

csr matrix_csr;
matrix_csr.nr = 3;
matrix_csr.nc = 5;
matrix_csr.nz = 8;

matrix_csr.row = malloc((matrix_csr.nr+1 * sizeof(int))); //Keeps track of where each rows starts ina flattened array of non zero elements of the matrix. we build it by checking total number of non zero entries see so far while processing rows. Refer to chatgpt chat for details on this discussion: https://chatgpt.com/share/69dcb3f7-5bdc-838d-aa63-8fe28d5b8432
// start with 0
// after row 0 → 2
// after row 1 → 5
// after row 2 → 8
// A.row = [0, 2, 5, 8]
matrix_csr.col = malloc((matrix_csr.nz * sizeof(int)));
matrix_csr.coeff = malloc(matrix_csr.nz * sizeof(double));
matrix_csr.coeff[0] = 3;
matrix_csr.coeff[1] = 1;
matrix_csr.coeff[2] = 1;
matrix_csr.coeff[3] = 4;
matrix_csr.coeff[4] = 5;
matrix_csr.coeff[5] = 3;
matrix_csr.coeff[6] = 2;
matrix_csr.coeff[7] = 4;

matrix_csr.col[0] = 1;
matrix_csr.col[1] = 2;
matrix_csr.col[2] = 0;
matrix_csr.col[3]= 2;
matrix_csr.col[4]=4;
matrix_csr.col[5]=1;
matrix_csr.col[6]=2;
matrix_csr.col[7]=3;

matrix_csr.row[0]=0;
matrix_csr.row[1]= 2;
matrix_csr.row[2]= 5;
matrix_csr.row[3] = 8;
double matrix_recont;


printf("\n");
//let's reconstruct row by  row. Just reconstruct first row for now.







double alpha = 2;
double beta = 3;
double vector_x[5] = {1,3,5,3,3};
double vector_y[5] = {1,4,2,4,2};
//problem b: CSR format.
matrix_mult_csr();
for(int i=0; i<matrix_csr.nr; i++){
	for(int j = matrix_csr.row[i]; j<matrix_csr.row[i+1]; j++){
		printf("%0.1f ", matrix_csr.coeff[j]);
		vector_y[i] += alpha*(matrix_csr.coeff[j])*(vector_x[matrix_csr.col[j]]);

}
printf("\n");
}
printf("\n");
for(int i=0; i<matrix_csr.nc; i++){
printf("vector_y[%d]: %f", i, vector_y[i]);
}
    return -1;
}