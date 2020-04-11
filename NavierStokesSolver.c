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
int itmax=30;
double h= 1/ (float)n;
//beta=omega*h^2/(4*dt);
double beta = 0.0542; 
double u[n+2][n+2] = {0};
double v[n+2][n+2] = {0}; 
double p[n+2][n+2] = {0};

int main()
{
    printf("%d, %d, %d, %f, %d, %f\n", n, Re, tmax, dt, itmax, h);  // just to test
    printf("%f\n", epsi);
    //int row, columns;
    for (int row=0; row<n+2; row++)
    {
        for(int columns=0; columns<n+2; columns++)
            {
            printf("%d     ", u[row][columns]);
            }
        printf("\n");
    }

    if (dt> fmin(fmin(h,Re*pow(h,2)/4),2/(float)Re)) {
        printf("Warning! dt should be less than %f ", fmin(fmin(h,Re*pow(h,2)/4),2/(float)Re));
    
    };
    double fux, fuy, fvx, fvy, visu, visv;
    
    for (float t=0.0; t<tmax; t= t+dt)
    {
        printf("%f\n", t);
        
        for(int i=2; i<n+1; i++)
        {
            for(int j=2; j<n+1; j++)
            {
                // update temporal veloceties for the internal nodes from prev timestep
                fux = (pow(u[i][j]+u[i+1][j],2) - pow(u[i-1][j]+u[i][j],2))*0.25/h;
                fuy = ((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1]) - (v[i][j-1]+v[i+1][j-1])*(u[i][j-1]+u[i][j]))*0.25/h;
                fvx = ((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j]) - (u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[i][j]))*0.25/h;
                fvy = (pow(v[i][j]+v[i][j+i],2) - pow(v[i][j-1]+v[i][j],2))*0.25/h;
                visu=(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-4.0*u[i][j])/(Re*pow(h,2));
                visv=(v[i+1][j]+v[i-1][j]+v[i][j+1]+v[i][j-1]-4.0*v[i][j])/(Re*pow(h,2));
                u[i][j]=u[i][j]+dt*((p[i][j] - p[i+1][j])/h - fux - fuy + visu);
                v[i][j]=v[i][j]+dt*((p[i][j] - p[i][j+1])/h - fvx - fvy + visv);
            }
        }
        
        for(int iter=1;iter<itmax; iter++)           // BcVel, Boundary conditions for the velocities
        {
            for(int j=1; j<n+2; j++)
            {
                //left wall
                u[1][j]=0.0;
                v[1][j]=-v[2][j];
                //right wall
                u[n+1][j]=0.0;
                v[n+2][j]=-v[n+1][j];
            }
            for(int i=1; i<n+2; i++)
            {
                //top wall (moving)
                u[i][n+2]=-u[i][n+1]+2.0;
                //printf("%f, ", u[i][n+2]);
                v[i][n+1]=0.0;
                //bottom wall
                u[i][1]=-u[i][2];
                v[i][1]=0.0;
            }
            //int iflag=0;               // Piter, Pressure iterations
            for (int j=2; j<n+1; j++)
            {
                for(int i=2; i<n+1; i++)
                {
                    double div=(u[i][j]-u[i-1][j])/h+(v[i][j]-v[i][j-1])/h;
                    //if (abs(div)>=epsi),iflag=1;end
                    double delp = -beta*div;
                    p[i][j]  =p[i][j]  +delp;
                    u[i][j]  =u[i][j]  +delp*dt/h;
                    u[i-1][j]=u[i-1][j]-delp*dt/h;
                    v[i][j]  =v[i][j]  +delp*dt/h;
                    v[i][j-1]=v[i][j-1]-delp*dt/h;
                }
            }
            //if(iflag==0)break,end
        }
    }
    for (int row=0; row<n+2; row++)
    {
        for(int columns=0; columns<n+2; columns++)
            {
            printf("%f     ", u[row][columns]);
            }
        printf("\n");
    }


    return 0;
}