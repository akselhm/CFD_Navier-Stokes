#include "header.h"

void cavityBc(double u[n+2][n+2], double v[n+2][n+2]){
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
        u[i][n+1] = -u[i][n] + 2;
        v[i][n]   = 0.0;
        //bottom wall
        u[i][0] = -u[i][1];
        v[i][0] = 0.0;
    }
}

void pipelineBc(double u[n+2][n+2], double v[n+2][n+2]){
    for(int j=0; j<n+2; j++)
    {
        //left side (inlet condition)
        u[0][j] = 1.0;
        v[0][j] = 0.0;
        //right side (outlet)
        u[n+1][j] = u[n][j];
    }
    for(int i=0; i<n+2; i++)
    {
        //top wall (no slip)
        u[i][n+1] = -u[i][n];
        v[i][n]   = 0.0;
        //bottom wall (no slip)
        u[i][0] = -u[i][1];
        v[i][0] = 0.0;
    }
}


void plateBc(double u[n+2][n+2], double v[n+2][n+2]){
    int platesize = n/4;

    int y0 = (n-platesize)/2;
    int y1 = (n+platesize)/2;

    for(int j=0; j<n+2; j++)
    {
        //left side (inlet condition)
        u[0][j] = 1.0;
        v[0][j] = 0.0;
        //right side (outlet)
        u[n+1][j] = u[n][j];
        v[n+1][j] = v[n][j]; 
    }
    for(int i=0; i<n+2; i++)
    {
        //top wall (free slip)
        v[i][n] = 0.0;
        u[i][n+1] = u[i][n];   
        //bottom wall (free slip)
        v[i][0] = 0.0;
        u[i][0] = u[i][1];  
    }
    //plate
    for(int j=y0; j<y1+1; j++)
    {
        //no slip condition on plate
        u[0][j] = 0.0;
        v[0][j] = -v[1][j];
    }
}

void boxBC(double u[n+2][n+2], double v[n+2][n+2], int boxLength, int boxHeigth, int boxStart){
    // x0, x1, y0, y1 are "coordinates" represented as int (node) between 0 og n so that they fit the closest grid cell
    //visualization of box(square/rectangle):
    //    (x0,y1)__ ______________(x1,y1)
    //          |__|           |__|
    //             |              |
    //    (x0,y0)__|______________|(x1,y0)
    //          |__|           |__|


    int x0 =   boxStart;
    int x1 = boxStart + boxLength;
    int y0 = (n - boxHeigth)/2;
    int y1 = (n + boxHeigth)/2;

    for(int j=0; j<n+2; j++)
    {
        //left side (inlet condition)
        u[0][j] = 1.0;
        v[0][j] = 0.0;
        //right side (outlet)
        u[n+1][j] = u[n][j];
        v[n+1][j] = v[n][j]; 
    }
    for(int i=0; i<n+2; i++)
    {
        //top wall (free slip)
        v[i][n] = 0.0;
        u[i][n+1] = u[i][n];   
        //bottom wall (free slip)
        v[i][0] = 0.0;
        u[i][0] = u[i][1]; 
    }
    //box
    for(int i= x0+1; i<=x0; i++)
    {
        for(int j=y0+1; j<=y1; j++)
        {
            // set all velocities inside and on top and backside to zero   
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
    for(int i= x0+1; i<x1; i++){
        // u = 0 along top and bottom (interpolated)
        u[i][y0+1] = -u[i][y0];
        u[i][y1] = -u[i][y1+1];
        // no flow in from below
        v[i][y0] = 0.0;

    }
    for(int j= y0+1; j<y1; j++){
        // v = 0 along front and backside (interpolated)
        v[x0+1][j] = -v[x0][j];
        v[x1][j] = -v[x1+1][j];
        // no flow in from the front
        u[x0][j] = 0.0;

    }
    u[x0][y1] = 0.0;
    v[x1][y0] = 0.0;
}