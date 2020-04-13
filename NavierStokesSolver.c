#include <stdio.h>
#include <stdlib.h>
#include <math.h>  

const int n = 10;             // must use const to use in u, v and p matrises
double epsi = pow(10,-6);     // pow(a,b) returnes a to the power of b
int nn[] = {0, 5, 10, 20, 30, 40, 60, 100, 500};
float oo[] = {1.7, 1.78, 1.86, 1.92, 1.95, 1.96, 1.97, 1.98, 1.99};
int Re = 100;
int tmax = 1;
double dt = 0.01;
int itmax=3;
double h= 1/ (float)n;
//beta=omega*h^2/(4*dt);
double beta = 0.0542; 
double u[n+2][n+2] = {0};
double v[n+2][n+2] = {0}; 
double p[n+2][n+2] = {0};


void updateInternal(double u[n+2][n+2], double v[n+2][n+2], int i, int j);  //function for updating temporal velocities on an internal node[i][j]
void updateBoundaries(double u[n+2][n+2], double v[n+2][n+2], int n);       //function for applying all boundary conditions
void pressureIt(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], int n, int iflag);  //function for pressure iterations and new velocities
void printDoubleMatrix(double matrix[n+2][n+2]);                            //function for printing a (n+2)x(n+2)-matrix

int main()
{
    double testmatrix[5][5] = {{1, 2.3, 4.2, 2.0, 3.6},
                        {2.4, 3.7, 1.1, 6, 3.2},
                        {0.2, 0.2, 0.2, 0.2, 0.2},
                        {3.5, 4.6, 2.2, 3, 1.1},
                        {2.4, 3.7, 1.1, 6, 3.2}};
    /*
    printf("%d, %d, %d, %f, %d, %f\n", n, Re, tmax, dt, itmax, h);  // just to test
    printf("%f\n", epsi);
    */
    if (dt> fmin(fmin(h,Re*pow(h,2)/4),2/(float)Re)) {
        printf("Warning! dt should be less than %f ", fmin(fmin(h,Re*pow(h,2)/4),2/(float)Re));
    
    };
    //double fux, fuy, fvx, fvy, visu, visv;
    
    for (float t=0.0; t<tmax; t= t+dt)
    {
        printf("%f\n", t);
        
        for(int i=1; i<n+1; i++)
        {
            for(int j=1; j<n+1; j++)
            {
                updateInternal(u, v, i, j);
            }
        }
        
        for(int iter=1;iter<itmax; iter++)           // BcVel, Boundary conditions for the velocities
        {
            updateBoundaries(u, v, n);
            int iflag=0;  
            pressureIt(u, v, p, n, iflag);
            if(iflag==0) 
            {
                break;
            }
        }
    }

    printf("u-matrix:\n");
    printDoubleMatrix(u);
    printf("\n");

    printf("v-matrix:\n");
    printDoubleMatrix(v);
    printf("\n");

    printf("p-matrix:\n");
    printDoubleMatrix(p);
    printf("\n");
    /*
    for(int i=1; i<4; i++)
    {
        for(int j=1; j<4; j++)
        {
            updateInternal(testmatrix, testmatrix, i,j);
        }
    }
    for (int row=0; row<5; row++)
    {
        for(int columns=0; columns<5; columns++)
            {
            printf("%f     ", testmatrix[row][columns]);
            }
        printf("\n");
    }*/
    return 0;
}

//must update to [n+2][n+2]
void updateInternal(double u[n+2][n+2], double v[n+2][n+2], int i, int j){
    // update temporal veloceties for the internal nodes from prev timestep
    double fux, fuy, fvx, fvy, visu, visv;
    fux = (pow(u[i][j]+u[i+1][j],2) - pow(u[i-1][j]+u[i][j],2))*0.25/h;
    fuy = ((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1]) - (v[i][j-1]+v[i+1][j-1])*(u[i][j-1]+u[i][j]))*0.25/h;
    fvx = ((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j]) - (u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[i][j]))*0.25/h;
    fvy = (pow(v[i][j]+v[i][j+i],2) - pow(v[i][j-1]+v[i][j],2))*0.25/h;
    visu=(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-4.0*u[i][j])/(Re*pow(h,2));
    visv=(v[i+1][j]+v[i-1][j]+v[i][j+1]+v[i][j-1]-4.0*v[i][j])/(Re*pow(h,2));
    u[i][j]=u[i][j]+dt*((p[i][j] - p[i+1][j])/h - fux - fuy + visu);
    v[i][j]=v[i][j]+dt*((p[i][j] - p[i][j+1])/h - fvx - fvy + visv);
}

void updateBoundaries(double u[n+2][n+2], double v[n+2][n+2], int n){
     for(int j=0; j<n+2; j++)
    {
        //left wall
        u[0][j]=0.0;
        v[0][j]=-v[1][j];
        //right wall
        u[n][j]=0.0;
        v[n+1][j]=-v[n][j];
    }
    for(int i=0; i<n+2; i++)
    {
        //top wall (moving)
        u[i][n+1]=-u[i][n]+2.0;
        v[i][n]=0.0;
        //bottom wall
        u[i][0]=-u[i][2];
        v[i][0]=0.0;
    }
}

void pressureIt(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], int n, int iflag) {             // Piter, Pressure iterations
    for (int j=1; j<n+1; j++)
    {
        for(int i=1; i<n+1; i++)
        {
            double div=(u[i][j]-u[i-1][j])/h+(v[i][j]-v[i][j-1])/h;
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