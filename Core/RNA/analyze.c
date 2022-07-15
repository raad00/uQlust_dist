//#include "stdafx.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "nrutil.h"
#include "rna.h"

double torsion(double **d)
/* get torsion angle a-b-c-d in degrees */
{
    double ang_deg, dij;
    double **vec3;
    long i, j;

    vec3 = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= 3; i++) {
        for (j = 1; j <= 3; j++)
            if (i == 1)
                vec3[i][j] = d[i][j] - d[i + 1][j];        /* b-->a */
            else
                vec3[i][j] = d[i + 1][j] - d[i][j];
        dij = veclen(vec3[i]);
        if (dij > BOND_UPPER_LIMIT) {
            free_dmatrix(vec3, 1, 3, 1, 3);
            return EMPTY_NUMBER;
        }
    }

    ang_deg = vec_ang(vec3[1], vec3[3], vec3[2]);
    free_dmatrix(vec3, 1, 3, 1, 3);

    return ang_deg;
}

