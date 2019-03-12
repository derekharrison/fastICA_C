/*
 * Setup.c
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#include "../inc/Memory.h"
#include "../inc/Parameter.h"
#include <stdio.h>


void setupVars()
{
    /*
     * This function sets up various parameters used by the algorithm
     */

    periodSource5 = (finalTime - initialTime)/na;
    periodSource6 = (finalTime - initialTime)/ns/2;

    timeVector     = matrix1D(M);

    Amix         = matrix2D(N, N);
    W             = matrix2D(N, N);
    WT             = matrix2D(N, N);

    /*Generating Data for ICA*/

    S             = matrix2D(N, M);
    Sest         = matrix2D(N, M);
    Xobs         = matrix2D(N, M);
    X             = matrix2D(N, M);
    Z             = matrix2D(N, M);
}
