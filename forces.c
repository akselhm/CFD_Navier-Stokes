#include "header.h"

void calforces(double p[n+2][n+2], int boxHeigth, int boxLength, int boxStart, double h){
    double lift = 0.0;
    double drag = 0.0;

    double p_top    = 0.0;
    double p_bottom = 0.0;
    double p_front  = 0.0;
    double p_behind = 0.0;

    int x0 = boxStart;
    int x1 = boxStart + boxLength;
    int y0 = (n - boxHeigth)/2;
    int y1 = (n + boxHeigth)/2;

    for(int i=x0+1; i<=x1; i++){
        lift += (p[i][y0] - p[i][y1+1])*h;
    }
    for(int j=y0+1; j<=y1; j++){
        drag += (p[x0][j] - p[x1+1][j])*h;
    }

    printf("Forces:\n");
    printf("lift:    %f\n", lift);
    printf("drag:    %f\n", drag);
}

