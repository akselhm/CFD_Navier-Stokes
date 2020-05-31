/*#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
//#include "skrivFunk.h"
#include <math.h> 
#include <time.h> 
*/
#include "header.h"
//#define n 14
#include <stdbool.h>



/*     // ------- lagt inn i header (gjort endringer, så ikke opdatert under) --------------
void updateInternal(double u[imax+2][jmax+2], double v[imax+2][jmax+2], int i, int j);  //function for updating temporal velocities on an internal node[i][j]
void updateBoundaries(double u[imax+2][jmax+2], double v[imax+2][jmax+2]);       //function for applying all boundary conditions
void pressureIt(double u[imax+2][jmax+2], double v[imax+2][jmax+2], double p[imax+2][jmax+2], int iflag);  //function for pressure iterations and new velocities

void printDoubleMatrix(double matrix[imax+2][jmax+2]);                            //function for printing (n+2)x(n+2)-matrix
void initializeMatrix(double matrix[imax+2][jmax+2], double initval);                // function for initializing a matrix such that all elements have value initval
void copyMatrix(double old[imax+2][jmax+2], double newM[imax+2][jmax+2]);               //update old to be identical to new

void makeBox(int istart, int iend, int jstart, int jend);
void boxBoundaries(double u[imax+2][jmax+2], double v[imax+2][jmax+2]);
*/

int main()
{
    double epsi = pow(10,-6);     // pow(a,b) returnes a to the power of b
    int nn[] = {0, 5, 10, 20, 30, 40, 60, 100, 500};
    float oo[] = {1.7, 1.78, 1.86, 1.92, 1.95, 1.96, 1.97, 1.98, 1.99};
    double Re = 100.0;
    int tmax = 10;
    double dt = 0.001;
    int itmax = 30;
    double h = 1.0/ (float)n;
    //beta=omega*h^2/(4*dt);  //omega finnes ved interpolasjon (ikke implementert enda)
    double beta = 0.0542; 
    double U_wall = 1.0;
    double u[n+2][n+2] = {0};
    double v[n+2][n+2] = {0}; 
    double p[n+2][n+2];

    double U[n][n];
    double V[n][n];
    double P[n][n];
    double psi[n+1][n+1] = {0};

    int boxplaced = 0;
/*
    struct Box{
        // a struct for a square/reqtangle-object
        // corners have coordinates: downleft=(istart, jstart), downright=(iend, jstart), upperright=(iend, jend), upperleft(istart, jend)
        int istart;
        int jstart;
        int iend;
        int jend;
    };
    struct Box box;
    //set coordinates of box
    box.istart = 3;
    box.iend = 7;
    box.jstart = 4;
    box.jend = 7;
*/

    //boxplaced = 1;
    /*
    printf("%d, %d, %d, %f, %d, %f\n", n, Re, tmax, dt, itmax, h);  // check if right variables
    printf("%f\n", epsi);
    */

    initializeMatrix(p, 1.0);

    if (dt> fmin(fmin(h,Re*pow(h,2)/4),2/(float)Re)) {
        printf("Warning! dt should be less than %f ", fmin(fmin(h,Re*pow(h,2)/4),2/(float)Re));
    
    };
    //double fux, fuy, fvx, fvy, visu, visv;
    if(boxplaced == 0){
        for (float t=0.0; t<tmax; t= t+dt)      //main loop
        {
            //printf("%f\n", t);
            
            for(int i=1; i<n+1; i++)
            {
                for(int j=1; j<n+1; j++)
                {
                    updateInternal(u, v, p, i, j, h, Re, dt);
                }
            }
            
            for(int iter=1;iter<itmax; iter++)           // BcVel, Boundary conditions for the velocities
            {
                updateBoundaries(u, v, U_wall);
                int iflag=0;  
                pressureIt(u, v, p, iflag, epsi, beta, h, dt);
                if(iflag==0) 
                {
                    break;
                }
            }
        }
    }
    /*
    else {      // ---square/reqtangle inside cavity ----
        for (float t=0.0; t<tmax; t= t+dt)      //main loop
        {
            //printf("%f\n", t);
            
            for(int i=1; i<n+1; i++)
            {
                for(int j=1; j<n+1; j++)
                {
                    if (i>=box.istart && i<=box.iend && j>=box.jstart && j<=box.jend) //inside box velocoties are 0 (apply no slip conditions later)
                    {
                        u[i][j] = 0.0;
                        v[i][j] = 0.0;
                    }
                    else if(i == box.istart-1 && j>=box.jstart && j<=box.jend){     //u is 0 from leftside into the box
                        u[i][j] = 0.0;
                    }
                    else if(j == box.jstart-1 && i>=box.istart && i<=box.iend){     //v is 0 from bottom into the box
                        v[i][j] = 0.0;
                    }
                    else{
                        updateInternal(n, u, v, p, i, j, h, Re);
                    }
                }
            }
            //v[box.istart]
            
            for(int iter=1;iter<itmax; iter++)           // BcVel, Boundary conditions for the velocities
            {
                updateBoundaries(n, u, v, U_wall);
                //boxBoundaries(u, v);
                int iflag=0;  
                pressureIt(n, u, v, p, iflag, epsi, beta, h, dt);
                if(iflag==0) 
                {
                    //break;
                }
            }
        }
    }
    */            // ---square/reqtangle inside cavity stop ----

    calculateUVP(u, v, p, U, V, P);
    streamlines(u, v, psi, h);

    // ---------- make vtf-file ---------
    int teller = 1;
    FILE *fidGL;
    fidGL = fopen("GLview.vtf", "w+");

    skrivHeader(fidGL);
    skrivRes(fidGL, & teller,P, psi, U, V);
    skrivSteg(fidGL, & teller);

    fclose(fidGL);
    // -----------------------------------

    // print velocity-matrices and pressure-matrix. The matrices are flipped 90 degrees clockwise 
    printf("\nu-matrix:\n");
    printDoubleMatrix(u);

    printf("\nv-matrix:\n");
    printDoubleMatrix(v);

    printf("\np-matrix:\n");
    printDoubleMatrix(p);

    return 0;
}




// ----------------   functions defined  -----------------------------------------------


/*
void updateInternal(int n, double u[n+2][n+2], double v[n+2][n+2], int i, int j){
    // update temporal veloceties for the internal nodes from prev timestep
    double fux, fuy, fvx, fvy, visu, visv;
    fux = ((u[i][j]+u[i+1][j])*(u[i][j]+u[i+1][j]) - (u[i-1][j]+u[i][j])*(u[i-1][j]+u[i][j]))*0.25/h;
    fuy = ((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1]) - (v[i][j-1]+v[i+1][j-1])*(u[i][j-1]+u[i][j]))*0.25/h;
    fvx = ((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j]) - (u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[i][j]))*0.25/h;
    fvy = ((v[i][j]+v[i][j+i])*(v[i][j]+v[i][j+i]) - (v[i][j-1]+v[i][j])*(v[i][j-1]+v[i][j]))*0.25/h;

    visu = (u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-4.0*u[i][j])/(Re*h*h);
    visv = (v[i+1][j]+v[i-1][j]+v[i][j+1]+v[i][j-1]-4.0*v[i][j])/(Re*h*h);
    
    u[i][j] = u[i][j]+dt*((p[i][j] - p[i+1][j])/h - fux - fuy + visu);
    v[i][j] = v[i][j]+dt*((p[i][j] - p[i][j+1])/h - fvx - fvy + visv);
}

void updateBoundaries(int n, double u[n+2][n+2], double v[n+2][n+2]){
    for(int j=0; j<n+2; j++)
    {
        //left wall
        u[0][j] = 0.0;
        v[0][j] = -v[1][j];
        //right wall
        u[n][j] = 0.0;
        v[n+1][j] = -v[n][j];
    }
    for(int i=0; i<n+2; i++)
    {
        //top wall (moving with the pace 1)
        u[i][n+1] = -u[i][n] + 2*U;
        v[i][n]   = 0.0;
        //bottom wall
        u[i][0] = -u[i][1];
        v[i][0] = 0.0;
    }
}

void pressureIt(int n, double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], int iflag) {             // Piter, Pressure iterations
    for (int j=1; j<n+1; j++)
    {
        for(int i=1; i<n+1; i++)
        {
            double div=(u[i][j]-u[i-1][j])/h+(v[i][j]-v[i][j-1])/h;     // continuity equation discretized with backwards differences
            if (fabs(div)>=epsi) {
                iflag=1;
            }
            double delp = -beta*div;
            p[i][j]  =p[i][j]  +delp;

            u[i][j]  =u[i][j]  +delp*dt/h;
            u[i-1][j]=u[i-1][j]-delp*dt/h;
            v[i][j]  =v[i][j]  +delp*dt/h;
            v[i][j-1]=v[i][j-1]-delp*dt/h;
        }
    }
}

*/
// ---------------------------------- UVP -----------------------------------------------------
/*
void calculateUVP(int n, double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], double U[n][n], double V[n][n], double P[n][n]){
    //flip the matrix 90 degrees anti clockwise and calculate velocities in the middle of the grid cells
    for (int i=0; i<n-1; i++)
    {
        for(int j=0; j<n-1; j++)
        {
            U[j][i] = (u[i][j+1]+u[i+1][i+1])/2;
            V[j][i] = (v[i+1][j]+v[i+1][i+1])/2;
            P[j][i] = p[i+1][j+1];
        }
    }
}
*/
// -----------------------------------streamlines ----------------------------------------------------------
/*
void streamlines(int n, double u[n+2][n+2], double v[n+2][n+2], double psi[n+1][n+1]){
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

    for (int i=0; i<n+1; i++)
    {
        for(int j=1; j<n+1; j++)
        {
            psi[j][i] = psi[i][j];
        }
    }
}
*/
// -------------------------matrix help-functions (print, copy, initialize) --------------------------------------
/*

void printDoubleMatrix(int n, double matrix[n+2][n+2]){
    for (int row=0; row<n+2; row++)
    {
        for(int columns=0; columns<n+2; columns++)
            {
            printf("%f     ", matrix[row][columns]);
            }
        printf("\n");
    }
}

void initializeMatrix(int n, double matrix[n+2][n+2], double initval) {
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

// ------------------------------------ Box ----------------------------------------------------
/*
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
        boxplaced = true;
    }
    
}

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