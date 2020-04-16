#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <time.h> 
#define n 10
#define imax 10
#define jmax 10

//const int n = 10;             // must use const to use in u, v and p matrises (can also use #define n 10)
double epsi = pow(10,-6);     // pow(a,b) returnes a to the power of b
int nn[] = {0, 5, 10, 20, 30, 40, 60, 100, 500};
float oo[] = {1.7, 1.78, 1.86, 1.92, 1.95, 1.96, 1.97, 1.98, 1.99};
float Re = 100.0;
int tmax = 10;
double dt = 0.001;
int itmax = 30;
double h = 1/ (float)imax;
//beta=omega*h^2/(4*dt);
double beta = 0.0542; 
double u[imax+2][jmax+2] = {0};
double v[imax+2][jmax+2] = {0}; 
double p[imax+2][jmax+2];


void updateInternal(double u[imax+2][jmax+2], double v[imax+2][jmax+2], int i, int j);  //function for updating temporal velocities on an internal node[i][j]
void updateBoundaries(double u[imax+2][jmax+2], double v[imax+2][jmax+2]);       //function for applying all boundary conditions
void pressureIt(double u[imax+2][jmax+2], double v[imax+2][jmax+2], double p[imax+2][jmax+2], int iflag);  //function for pressure iterations and new velocities
void printDoubleMatrix(double matrix[imax+2][jmax+2]);                            //function for printing (n+2)x(n+2)-matrix
void initializeMatrix(double matrix[imax+2][jmax+2], double initval);                // function for initializing a matrix such that all elements have value initval

int main()
{
    /*
    printf("%d, %d, %d, %f, %d, %f\n", n, Re, tmax, dt, itmax, h);  // just to test
    printf("%f\n", epsi);
    */

    initializeMatrix(p, 1.0);

    if (dt> fmin(fmin(h,Re*pow(h,2)/4),2/(float)Re)) {
        printf("Warning! dt should be less than %f ", fmin(fmin(h,Re*pow(h,2)/4),2/(float)Re));
    
    };
    //double fux, fuy, fvx, fvy, visu, visv;
    
    for (float t=0.0; t<tmax; t= t+dt)      //main loop
    {
        //printf("%f\n", t);
        
        for(int i=1; i<imax+1; i++)
        {
            for(int j=1; j<jmax+1; j++)
            {
                updateInternal(u, v, i, j);
            }
        }
        
        for(int iter=1;iter<itmax; iter++)           // BcVel, Boundary conditions for the velocities
        {
            updateBoundaries(u, v);
            int iflag=0;  
            pressureIt(u, v, p, iflag);
            if(iflag==0) 
            {
                break;
            }
        }
    }

    // print velocity-matrices and pressure-matrix. The matrices is flipped 90 degrees clockwise 
    printf("\nu-matrix:\n");
    printDoubleMatrix(u);

    printf("\nv-matrix:\n");
    printDoubleMatrix(v);

    printf("\np-matrix:\n");
    printDoubleMatrix(p);
    /*
    // OUTPUT DATA
	FILE *fout2, *fout3;
	fout2 = fopen("UVP.plt","w+t");
	fout3 = fopen("Central_U.plt","w+t");

	if ( fout2 == NULL )
	{
        printf("\nERROR when opening file\n");
        fclose( fout2 );
	}

    else
	{
	    fprintf( fout2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
	    fprintf( fout2, "ZONE  F=POINT\n");
	    fprintf( fout2, "I=%d, J=%d\n", imax+2, jmax+2 );

        for (int j = 0 ; j < (jmax+2) ; j++ )
        {
            for (int i = 0 ; i < (imax+2) ; i++ )
            {
                double xpos, ypos;
                xpos = i*h;
                ypos = j*h;

                fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, u[i][j], v[i][j], p[i][j] );
            }
        }
	}

	fclose( fout2 );
	
	// CENTRAL --U
    fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
    fprintf(fout3, "ZONE F=POINT\n");
    fprintf(fout3, "I=%d\n", n+2 );

    for (int j = 0 ; j < n+2 ; j++ )
    {
        double ypos;
        ypos = (double) j*h;

        fprintf( fout3, "%5.8lf\t%5.8lf\n", (u[(n+2)/2][j] + u[((n+2)/2)+1][j])/(2.), ypos );
    }
    fclose(fout3);
*/
    return 0;
}

void updateInternal(double u[imax+2][jmax+2], double v[imax+2][jmax+2], int i, int j){
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

void updateBoundaries(double u[imax+2][jmax+2], double v[imax+2][jmax+2]){
     for(int j=0; j<jmax+2; j++)
    {
        //left wall
        u[0][j] = 0.0;
        v[0][j] = -v[1][j];
        //right wall
        u[imax][j] = 0.0;
        v[imax+1][j] = -v[n][j];
    }
    for(int i=0; i<imax+2; i++)
    {
        //top wall (moving with the pace 1)
        u[i][jmax+1] = -u[i][jmax]+2.0;
        v[i][jmax]   = 0.0;
        //bottom wall
        u[i][0] = -u[i][1];
        v[i][0] = 0.0;
    }
}

void pressureIt(double u[imax+2][jmax+2], double v[imax+2][jmax+2], double p[imax+2][jmax+2], int iflag) {             // Piter, Pressure iterations
    for (int j=1; j<jmax+1; j++)
    {
        for(int i=1; i<imax+1; i++)
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

void printDoubleMatrix(double matrix[imax+2][jmax+2]){
    for (int row=0; row<imax+2; row++)
    {
        for(int columns=0; columns<jmax+2; columns++)
            {
            printf("%f     ", matrix[row][columns]);
            }
        printf("\n");
    }
}

void initializeMatrix(double matrix[imax+2][jmax+2], double initval) {
    for (int row=0; row<imax+2; row++)
    {
        for(int columns=0; columns<jmax+2; columns++)
        {
            matrix[row][columns] = initval;
        }
    }
}