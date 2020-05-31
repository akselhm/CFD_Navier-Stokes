#include "header.h"

void updateInternal(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], int i, int j,double h, double Re, double dt){
    // update temporal veloceties for the internal nodes from prev timestep and prev pressure
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


void updateBoundaries(double u[n+2][n+2], double v[n+2][n+2], double U_wall){
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
        //top wall (moving with the pace U_wall)
        u[i][n+1] = -u[i][n] + 2*U_wall;
        v[i][n]   = 0.0;
        //bottom wall
        u[i][0] = -u[i][1];
        v[i][0] = 0.0;
    }
}


void pressureIt(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], int iflag, double epsi, double beta, double h, double dt) {             // Piter, Pressure iterations
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