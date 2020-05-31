#include "header.h"

void skrivHeader(FILE *f){
    fprintf(f, "*VTF-1.00 \n!noder \n*NODES         1"); //Skriver kordinate til nodene
    for (int i=0; i<n; i=i+1){
        for( int j=0; j<n; j=j+1){
            fprintf(f, "\n%f   %f   %f", (j+1.0), (i+1.0), 0.0);
        }
    }
    
    fprintf(f, "\n\n*ELEMENTS            1\n%%NODES #1\n%%QUADS"); //Skriver nodene som former hvert element
    for (int j=0; j<(n-1); j=j+1){
        for (int i=1; i<n; i=i+1){
            fprintf(f, "\n%i   %5i    %5i    %5i", i+(j*n), i+1+(j*n), i+n+1+(j*n), i+n+(j*n));
        }
    }
    fprintf(f, "\n\n*GLVIEWGEOMETRY 1\n%%ELEMENTS\n1\n");
}

void skrivRes(FILE *f, int *teller, double P[n][n], double psi[n+1][n+1], double U[n][n],
              double V[n][n]){
    fprintf(f, "\n\n*RESULTS   %i\n%%DIMENSION 1\n%%PER_NODE #1",*teller); //Skriver trykket til fil
    for (int i=0; i<n; i=i+1){
        for( int j=0; j<n; j=j+1){
            fprintf(f, "\n%f", P[i][j]);
        }
    }
    
    *teller=*teller+1;
    
    fprintf(f, "\n\n*RESULTS   %i\n%%DIMENSION 1\n%%PER_NODE #1",*teller); //Skriver strømfunksjon til fil
    for (int i=0; i<n; i=i+1){
        for( int j=0; j<n; j=j+1){
            fprintf(f, "\n%f", psi[j][i]);
        }
    }
    
    *teller=*teller+1;
    
    fprintf(f, "\n\n*RESULTS   %i\n%%DIMENSION 3\n%%PER_NODE #1",*teller); //Skriver hastighetene til fil
    for (int i=0; i<n; i=i+1){
        for( int j=0; j<n; j=j+1){
            fprintf(f, "\n%f   %f   %f", U[i][j], V[i][j],0.0);
        }
    }
    
    *teller=*teller+1;
}

void skrivSteg(FILE *f, int *teller){
    fprintf(f,"\n\n*GLVIEWSCALAR 1\n%%NAME \"Pressure\" "); //Skriver trykkstegene
    for (int i=1; i<=(*teller/3); i=i+1){
        fprintf(f,"\n%%STEP  %i \n%i",i,(i*3)-2);
    }
    
    fprintf(f,"\n\n*GLVIEWSCALAR 2\n%%NAME \"Streamfunction\" "); //Skriver strømfunksjonstegene
    for (int i=1; i<=(*teller/3); i=i+1){
        fprintf(f,"\n%%STEP  %i \n%i",i,(i*3)-1);
    }
    
    fprintf(f,"\n\n*GLVIEWVECTOR 1\n%%NAME \"Velocity\" "); //Skriver hastighetsstegene
    for (int i=1; i<=(*teller/3); i=i+1){
        fprintf(f,"\n%%STEP  %i \n%i",i,i*3);
    }
}
