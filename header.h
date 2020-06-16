
#pragma once 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#define n 25

// ----------- loopfuncs.c ----------------
void updateInternal(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], int i, int j,double h, double Re, double dt); 
void pressureIt(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], int *iflag, double epsi, double beta, double h, double dt, int option, int boxHeigth, int boxStart, int boxLength); 

// ------------ BCs.c ---------------------
void cavityBc(double u[n+2][n+2], double v[n+2][n+2]); 
void pipelineBc(double u[n+2][n+2], double v[n+2][n+2]);
void plateBc(double u[n+2][n+2], double v[n+2][n+2]);
void boxBC(double u[n+2][n+2], double v[n+2][n+2], int boxLength, int boxHeigth, int boxStart);

// ------------ calculateUVP.c -------------
void calculateUVP(double u[n+2][n+2], double v[n+2][n+2], double p[n+2][n+2], double U[n][n], double V[n][n], double P[n][n]);

// ------------ streamlines.c ---------------
void streamlines(double u[n+2][n+2], double v[n+2][n+2], double psi[n+1][n+1], double h);

// ------------ forces.c --------------------
void calforces(double p[n+2][n+2], int boxHeigth, int boxLength, int boxStart, double h);

// ------------ matrixhelpfunctions.c -------- 
void printDoubleMatrix(double matrix[n+2][n+2]);
void checkforNaN(double matrix[n+2][n+2], int *hasNaN);

// ------------ visualization.c ---------------
void skrivHeader(FILE *f);
void skrivRes(FILE *f, int *teller, double P[n][n], double psi[n+1][n+1], double U[n][n], double V[n][n]);
void skrivSteg(FILE *f, int *teller);


