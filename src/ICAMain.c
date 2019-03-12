/*
 * ICAMain.c
 *
 *  Created on: 14 apr. 2016
 *      Author: Derek W. Harrison
 *
 *  This code is a C implementation of the FastICA method for recovering independent
 *  components in component mixtures
 *
 *  For more information refer to the paper 'ICA: Algorithms and applications' which is found
 *  in the same directory as the source files of this code.
 */

#include "../inc/../inc/Parameter.h"
#include "../inc/Functions.h"
#include "../inc/MatrixOps.h"
#include "../inc/AlgoFunctions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main()
{
    srand((unsigned)time(NULL));

    begin = clock();

    /*User input parameter data*/
    ParameterInput();

    /*Setting up variables and generating Data*/
    setupVars();

    /*Initializing mixing matrix*/
    Initialize();

    /*Setting up source signals*/
    SetUpSources();

    /*Generating observed sample data*/
    Xobs = XobsGen(Amix, S, N, M);

    /*FastICA algorithm*/
    X = PreprocessingCentering(Xobs, N, M);

    Z = PreprocessingWhitening(X, N, M);

    W = FastICA(Z, N, M, iterations);

    /*Outputting results of FastICA algorithm*/
    WT   = MatTranspose(W, N, N);

    Sest = MatMult(WT, N, N, Z, N, M);

    printMatrixS(N, M, Sest);

    /*Exporting estimated source data to .txt for visualization in Matlab*/
    ExportingData();

    /*Cleaning up*/
    FreeMemory();

    end         = clock();
    time_spent     = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("\ntimespent: %f\n", time_spent);

    return 0;
}
