#include "header.h"

void printDoubleMatrix(double matrix[n+2][n+2]){
    for (int row=0; row<n+2; row++)
    {
        for(int columns=0; columns<n+2; columns++)
            {
            printf("%f     ", matrix[row][columns]);
            }
        printf("\n");
    }
}


void checkforNaN(double matrix[n+2][n+2], int *hasNaN){
    for (int row=0; row<n+2; row++)
    {
        for(int columns=0; columns<n+2; columns++)
        {
            if (matrix[row][columns] != matrix[row][columns]){
                //printDoubleMatrix(matrix);
                *hasNaN = 1;
            }
        }
    }
}
/*  // ------------ Excessive helpfunctions ------------------

void initializeMatrix(double matrix[n+2][n+2], double initval) {
    for (int row=0; row<n+2; row++)
    {
        for(int columns=0; columns<n+2; columns++)
        {
            matrix[row][columns] = initval;
        }
    }
}

void copyMatrix(int n, double old[n+2][n+2], double newM[n+2][n+2]){
    for (int row=0; row<n+2; row++)
    {
        for(int columns=0; columns<n+2; columns++)
            {
            old[row][columns] = newM[row][columns];
            }
    }
}
*/