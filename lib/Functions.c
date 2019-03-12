/*
 * Functions.c
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../inc/Parameter.h"
#include "../inc/Functions.h"
#include "../inc/Memory.h"


void Initialize()
{
    /*
     * This function generates the mixing matrix Amix, used to
     * generate the observation matrix Xobs, containing the M observed samples of
     * N sources.
     */

    int i, j;
    /*Generating mixing matrix Amix*/
    for(i = 0; i < N; ++i)
        for(j = 0; j < N; ++j)
        {
            Amix[i][j] = (double)rand()/((double)RAND_MAX);
        }
}


void SetUpSources()
{
    /*
     * This function generated the source matrix S, containing the sources used
     * to generate data for the observation matrix.
     */

    int i;
    for(i = 0; i < M; ++i)
        timeVector[i] = (finalTime-initialTime)/(M-1)*i;

    /*Source 1*/
    for(i = 0; i < M; ++i)
        S[0][i] = funcSource1(timeVector[i]);

    /*Source 2*/
    for(i = 0; i < M; ++i)
        S[1][i] = funcSource2(timeVector[i]);

    /*Source 3*/
    for(i = 0; i < M; ++i)
        S[2][i] = funcSource3(timeVector[i]);

    /*Source 4*/
    for(i = 0; i < M; ++i)
        S[3][i] = funcSource4(timeVector[i]);

    /*Source 5*/
    avgsource5 = 0.0;
    for(i = 0; i < M; ++i)
    {
        S[4][i] = funcSource5(timeVector[i]);
        avgsource5 += S[4][i];
    }

    /*Averaging Source 5*/
    avgsource5/=M;
    for(i = 0; i < M; ++i)
        S[4][i] -= avgsource5;

    /*Source 6*/
    avgsource6 = 0.0;
    for(i = 0; i < M; ++i)
    {
        S[5][i] = funcSource6(timeVector[i]);
        avgsource6 += S[5][i];
    }

    /*Averaging Source 6*/
    avgsource6/=M;
    for(i = 0; i < M; ++i)
        S[5][i] -= avgsource6;
}


double** XobsGen(double** Amix, double** S, int N, int M)
{
    /*
     * This function generates the observation matrix Xobs
     * consisting of M observed mixture samples of size N.
     */

    int i, j, k;
    double** Xobs = matrix2D(N, M);

    /*Generating observation matrix Xobs*/
    for(i = 0; i < N; ++i)
        for(j = 0; j < M; ++j)
        {
            for(k = 0; k < N; ++k)
                Xobs[i][j] += Amix[i][k]*S[k][j];
        }

    return Xobs;
}


void printMatrixS(int N, int M, double **S)
{
    /*
     * This function prints out the contents of matrices transposed.
     */

    int i, j;
    for(i = 0; i < M; ++i)
    {
        for(j = 0; j < N; ++j)
        {
            printf("%f\t", S[j][i]);
        }
        printf("\n");
    }
}


void printMatrix(int N, int M, double **S)
{
    /*
     * This function prints out the contents of matrices.
     */

    int i, j;

    for(i = 0; i < N; ++i)
    {
        for(j = 0; j < M; ++j)
        {
            printf("%f\t", S[i][j]);
        }
        printf("\n");
    }
}

void FreeMemory()
{
    /*
     * This function frees memory allocated for execution of the
     * algorithm.
     */

    free(Amix); free(W); free(WT);
    free(timeVector);
    free(S);
    free(Sest);
    free(Xobs);
    free(X);
    free(Z);
}
