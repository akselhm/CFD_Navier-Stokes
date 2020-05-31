/*#include "header.h"


    struct Box{
        // a struct for a square/reqtangle-object
        // corners have coordinates: downleft=(istart, jstart), downright=(iend, jstart), upperright=(iend, jend), upperleft(istart, jend)
        int istart;
        int jstart;
        int iend;
        int jend;
    };

void makeBox(int istart, int iend, int jstart, int jend){
    if (istart >= iend || jstart >= jend){
        printf("Unvalid inputvalues");
    }
    else if (istart < 1 || jstart < 1 || jstart < 1 || jstart < 1){
        printf("Unvalid inputvalues");
    }
    else{
        struct Box box;
        box.istart = istart;
        box.iend =iend;
        box.jstart = jstart;
        box.jend = jend;
        //boxplaced = true;
    }
    
}
*/

/*
void boxBoundaries(int n, double u[n+2][n+2], double v[n+2][n+2]){
    int i0 = 3;//box->istart;
    int i1 = 7;//box->iend;
    int j0 = 4;//box->jstart;
    int j1 = 7;//box->jend;
    // velocities 0.0 inside the box
    for(int i=i0+1; i<i1; i++)
    {
        for(int j=i0+1; j<i1; i++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
    //boundaries
    for(int j=i0; j<i1+1; j++)
    {
        //left wall
        u[i0-1][j] = 0.0;
        v[i0][j] = -v[i0-1][j];
        //right wall
        u[i1][j] = 0.0;
        v[i1][j] = -v[i1+1][j];
    }
    for(int i=i0; i<i1+1; i++)
    {
        //top wall (moving with the pace 1)
        u[i][j1] = -u[i][j1+1];
        v[i][j1]   = 0.0;
        //bottom wall
        u[i][j0] = -u[i][j0-1];
        v[i][j0] = 0.0;
    }
}
*/