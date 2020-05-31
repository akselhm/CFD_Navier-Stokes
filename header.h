#pragma once 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#define n 14

// --------------------------------- NavierStokesSolver --------------------------------------------------
void updateInternal(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], int i, int j,double h, double Re, double dt);  //function for updating temporal velocities on an internal node[i][j]
void updateBoundaries(double u[n+2][n+2], double v[n+2][n+2], double U_wall);       //function for applying all boundary conditions
void pressureIt(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], int iflag, double epsi, double beta, double h, double dt);  //function for pressure iterations and new velocities

void calculateUVP(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], double U[n][n], double V[n][n], double P[n][n]);
void streamlines(double u[n+2][n+2], double v[n+2][n+2], double psi[n+1][n+1], double h);

void printDoubleMatrix(double matrix[n+2][n+2]);                            //function for printing (n+2)x(n+2)-matrix
void initializeMatrix(double matrix[n+2][n+2], double initval);                // function for initializing a matrix such that all elements have value initval
//void copyMatrix(int n, double old[n+2][n+2], double newM[n+2][n+2]);               //update old to be identical to new

//void makeBox(int istart, int iend, int jstart, int jend);
//void boxBoundaries(int n, double u[n+2][n+2], double v[n+2][n+2]);

// --------------------------------- skrivHeader --------------------------------------------------

void skrivHeader(FILE *f);

// --------------------------------- skrivRes --------------------------------------------------
void skrivRes(FILE *f, int *teller, double P[n][n], double psi[n+1][n+1], double U[n][n], double V[n][n]);

// --------------------------------- skrivSteg --------------------------------------------------
void skrivSteg(FILE *f, int *teller);


