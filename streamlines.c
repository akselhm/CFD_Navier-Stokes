#include "header.h"

void streamlines(double u[n+2][n+2], double v[n+2][n+2], double psi[n+1][n+1], double h){
    for(int i=1; i<n+1; i++)
    {
        psi[i][0] = - h*v[i][0] + psi[i-1][0];
    }

    for (int i=0; i<n+1; i++)
    {
        for(int j=1; j<n+1; j++)
        {
            psi[i][j] = h*u[i][j] + psi[i][j-1];
        }
    }
    // endre siste her til å rotere matrisen
    for (int i=0; i<n+1; i++)
    {
        for(int j=1; j<n+1; j++)
        {
            psi[j][i] = psi[i][j];
        }
    }
}

