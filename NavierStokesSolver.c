#include "header.h"


int main()
{
    double epsi = 1e-6; 
    int nn[] = {0, 5, 10, 20, 30, 40, 60, 100, 500};
    float oo[] = {1.7, 1.78, 1.86, 1.92, 1.95, 1.96, 1.97, 1.98, 1.99};
    double Re = 80.0;
    int tmax = 20;
    double dt = 0.0001;
    int itmax = 2000;
    double h = 1.0/ (float)n;
    double omega;   //estimate later
    double u[n+2][n+2] = {0.0};
    double v[n+2][n+2] = {0.0}; 
    double p[n+2][n+2] = {0.0};

    double U[n][n];
    double V[n][n];
    double P[n][n];
    double psi[n+1][n+1] = {0.0};

    int hasNaN = 0;     // update to 1 if a matrix has a NaN

    int option = 2; 
    /*
        Valid options:
        1 = Flow in cavity box
        2 = Flow around square
        3 = Flow in pipeline
        4 = Flow behind plate
    */    
    if (option < 1 || option > 4) {
        printf("Warning! Unvalid option. Option must be 1, 2, 3 or 4");
        return 0;
    }

    // Box properties --------------------
    int boxLength = n/8;
    int boxHeigth = n/8;
    int boxStart = n/3;  
    // --------------------------
    
    
    // --------------- estimate omega and find beta -----
    for (int i = 1; i<sizeof(nn); i++){
        if (n < nn[i]){
            omega = oo[i-1] + (oo[i] - oo[i-1]) * (n - nn[i-1])/(nn[i]-nn[i-1]);
            break;
        }
    }
    double beta = omega*h*h/(4*dt);
    
    printf("%f\n", omega);
    printf("%f\n", beta);
    // ------------------ check stability ---------------------------------
    double criticalvalue = fmin(fmin(h, Re*h*h/4), 2.0/Re);
    if (dt > criticalvalue) {
        printf("Warning! dt should be less than %f ", criticalvalue);
        return 0;
    }

    for (float t=0.0; t<tmax; t= t+dt)      //main loop
    {   
        for(int i=1; i<n+1; i++)
        {
            for(int j=1; j<n+1; j++)
            {
                if (option != 2 || j <= (n - boxHeigth)/2 || j > (n  + boxHeigth)/2 || i <= boxStart || i > boxStart + boxLength){
                    // update temporal velocities based on old pressre with the Navier Stokes equations if cell_i,j not inside square/rectangle
                    updateInternal(u, v, p, i, j, h, Re, dt);
                }
            }
        }
        
        for(int iter=1;iter<itmax; iter++)     
        {
            if(option == 1)
            {
                cavityBc(u, v);
            }
            else if(option == 2)
            {
                boxBC(u, v, boxLength, boxHeigth, boxStart);
            }
            else if(option == 3)
            {
                pipelineBc(u,v);
            }   
            else if(option == 4)
            {
                plateBc(u,v);
            }

            int iflag=0;  
            pressureIt(u, v, p, &iflag, epsi, beta, h, dt, option, boxHeigth, boxStart, boxLength);
            if(iflag==0) 
            {
                break;
            }
        }
        checkforNaN(p, &hasNaN); //checks in p-matrix, as it is considered to be most likely to first obtain a NaN
        if (hasNaN == 1){
            printf("Warning! One or more element in one of the matrices is NaN at time %f. Change parameters and try again.", t);
            return 0;
        }
    }

    calculateUVP(u, v, p, U, V, P);
    streamlines(u, v, psi, h);

    if(option == 2){
        calforces(p, boxHeigth, boxLength, boxStart, h);
    }
    // ---------- make vtf-file ---------
    int teller = 1;
    FILE *fidGL;
    fidGL = fopen("GLview.vtf", "w+");

    skrivHeader(fidGL);
    skrivRes(fidGL, & teller,P, psi, U, V);
    skrivSteg(fidGL, & teller);

    fclose(fidGL);
    // -----------------------------------

    // print velocity-matrices and pressure-matrix (for testing purposes). The matrices are flipped 90 degrees clockwise 
/*
    printf("\nu-matrix:\n");
    printDoubleMatrix(u);

    printf("\nv-matrix:\n");
    printDoubleMatrix(v);

    printf("\np-matrix:\n");
    printDoubleMatrix(p);
*/
    return 0;
}
