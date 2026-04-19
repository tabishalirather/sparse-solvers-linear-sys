#include <stdio.h>
#include <stdlib.h>
// #include "../../Library/csr_mv.h"

int main() {
    printf("What is happening?\n");
    typedef struct csr_matrix {
        int nr;
        int nc;
        int nz;
        int *row;
        int *col;
        double *coeff;
    } csr;

    csr A;

    // Basic info
    A.nr = 3;
    A.nc = 5;
    A.nz = 8;

    // Allocate memory
    A.row = (int*) malloc((A.nr + 1) * sizeof(int));
    A.col = (int*) malloc(A.nz * sizeof(int));
    A.coeff = (double*) malloc(A.nz * sizeof(double));

    // Fill row pointer array
    A.row[0] = 0;
    A.row[1] = 2;
    A.row[2] = 5;
    A.row[3] = 8;

    // Fill column indices
    A.col[0] = 1;
    A.col[1] = 2;
    A.col[2] = 0;
    A.col[3] = 2;
    A.col[4] = 4;
    A.col[5] = 1;
    A.col[6] = 2;
    A.col[7] = 3;

    // Fill coefficients
    A.coeff[0] = 3;
    A.coeff[1] = 1;
    A.coeff[2] = 1;
    A.coeff[3] = 4;
    A.coeff[4] = 5;
    A.coeff[5] = 3;
    A.coeff[6] = 2;
    A.coeff[7] = 4;

    // float a = 0.5;
    // float b = 2;
    // double vector_x[3] = {4,2,1};
    // double vector_y[5] = {2,6,1};
    // double y=0;

    // This is wrong because it mistmatches the indicecs, just removing zero does nto solve the problem.
    // for (int i=0; i<A.nr; i++) {
    //     int k = 0;
    //     for (int j=A.row[i]; j<A.row[i+1]; j++) {
    //
    //
    //         //The inner for loop runs like this:
    //         // first run: from 0 to 2 which is 0th and 1th iteration
    //         // second run: from 2 to 5 whch is 2nd, 3rd and 4th iteration, so on and so forth. So depending on the difference between n+1th and nth element in A.row it assings elements to the same row, until it skips to the next row. So in the first iteration, the first two elements are assigned to the first row, in the second iteration, then next three elements from A.coeff are assigned to second row, and given that we already have column index of non zero elements, we can reconstruct the array.
    //         //printf("\n(i,j): (%d, %d)\n ", i,j);
    //         printf("%f ", a*A.coeff[j]*vector_x[k] + b*vector_y[k]);
    //         //Implement y <- a*Ax + b*y
    //         y += a*A.coeff[j]*vector_x[k]+b*vector_y[k];
    //         k++;
    //        // printf("A is: %f \n", A.coeff[j]);
    //         //printf("x is: %f \n", vector_x[k++]);
    //     }
    //     printf("\n");
    // }
    // printf("y is: %f", y);

    float a = 1;
    float b = 1;
    double vector_x[5] = {1,7,-2,4,1}; //Same num of rows as the columns original matrix with which it is multiplied, (not the CSR reduced format).
    double vector_y[3] = {2,6,1}; // this should be same as the rows of original matrix,( not in the reduced format)
    /*Here A is 3*5, so X has to be 5*1 and Ax will be 3*1 so y has to be 3*1*/
    /*
    * // CSR matrix-vector multiplication:
    // For each row, iterate over its nonzero entries.
    // A.col[j] tells which x element to use.
    // Compute sum += A.coeff[j] * x[A.col[j]] (ignore zeros).
     */
    for (int i=0; i<A.nr; i++) {
        double sum = 0;
        for (int j=A.row[i]; j<A.row[i+1]; j++) {
            sum += A.coeff[j] * vector_x[A.col[j]];
        }
        vector_y[i] = a*sum + b*vector_y[i];
    }

    printf("Results of matrix mult are: \n");
    for (int i=0; i<A.nr; i++) {
        printf("%f\n", vector_y[i]);
    }
    free(A.col);
    free(A.row);
    free(A.coeff);
return -1;
}
