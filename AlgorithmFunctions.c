/*
 * AlgorithmFunctions.c
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#include "MatrixOps.h"
#include "Memory.h"
#include <math.h>
#include <stdlib.h>


double** PreprocessingCentering(double** Xobs, int N, int M)
{
	/*
	 * This function performs the centering operation on a matrix X.
	 * For more details refer to the reference literature (ICA: Algorithms and applications).
	 */

	int i, j;
	double **X, *meanVector;

	X			= matrix2D(N, M);
	meanVector 	= matrix1D(N);

	/*Calculating mean vector of observation matrix*/
	for(i = 0; i < N; ++i)
		for(j = 0; j < M; ++j)
			meanVector[i] += Xobs[i][j];

	for(i = 0; i < N; ++i)
		meanVector[i] /= M;

	/*Centering observation matrix Xobs*/
	for(i = 0; i < N; ++i)
		for(j = 0; j < M; ++j)
			X[i][j] = Xobs[i][j] - meanVector[i];

	free(meanVector);

	return X;
}


double** PreprocessingWhitening(double** X, int N, int M)
{
	/*
	 * This function performs the whitening operation on a matrix X.
	 * For more details refer to the reference literature (ICA: Algorithms and applications)
	 * Note: the number of iterations for used for eigen decomposition of matrix ExxT is by default 100
	 * The default value can be changed by changing iterationsED to desired value.
	 */

	int i, j, k;
	double  *EigValues, **EigVectors, **EigVectorsT, *Dnegroot;
	double **ExxT, **Z, **Dummy1, **Dummy2, **Drootmat;

	int iterationsED = 100;

	ExxT 		= matrix2D(N, N);
	EigValues 	= matrix1D(N);
	Dnegroot	= matrix1D(N);
	Drootmat 	= matrix2D(N, N);
	Dummy1 		= matrix2D(N, N);
	Dummy2 		= matrix2D(N, N);
	EigVectors	= matrix2D(N, N);
	EigVectorsT	= matrix2D(N, N);
	Z 			= matrix2D(N, M);

	/*Calculating covariance*/
	for(i = 0; i < N; ++i)
		for(j = 0; j < N; ++j)
		{
			for(k = 0; k < M; ++k)
				ExxT[i][j] += X[i][k]*X[j][k] / (M - 1);
		}

	/*Eigen Decomposition of (N x N real symmetric) covariance matrix ExxT of X*/
	EigenDecomposition(ExxT, N, EigVectors, EigValues, iterationsED);

	/*Building matrix D^-1/2, containing the negative roots of the eigenvalues of covariance matrix of X,
	 * the centered data of Xobs*/
	for(i = 0; i < N; ++i)
		Dnegroot[i] = 1/sqrt(EigValues[i]);

	for(i = 0; i < N; ++i)
		for(j = 0; j < N; ++j)
		{
			if(i == j)
				Drootmat[i][j] = Dnegroot[i];
			else
				Drootmat[i][j] = 0.0;
		}

	/*Whitening matrix X. Z = E*1/sqrt(D)*E'*X */
	/*Dummy1 = E*1/sqrt(D)*/
	Dummy1 = MatMult(EigVectors, N, N, Drootmat, N, N);

	/*Tranpose of E*/
	EigVectorsT = MatTranspose(EigVectors, N, N);

	/*Dummy2 = Dummy1*E'*/
	Dummy2 = MatMult(Dummy1, N, N, EigVectorsT, N, N);

	/*Whitened matrix Z*/
	Z = MatMult(Dummy2, N, N, X, N, M);

	/*Clean Memory, unsure if actually needed*/
	free(EigValues);
	free(Dnegroot);
	free(Dummy1); free(Dummy2);
	free(Drootmat);
	free(EigVectors);
	free(EigVectorsT);
	free(ExxT);

	return Z;
}


double** FastICA(double ** Z, int N, int M, int iterations)
{
	/*
	 * This function performs the FastICA algorithm.
	 * For more details refer to the reference literature (ICA: Algorithms and applications).
	 */

	int 	p, i, j, k, it;
	double 	*G, *Gder, *ZGt, *dumsum, **W, *wp;;
	double 	GderOnes, f;

	G 			= matrix1D(M);
	Gder 		= matrix1D(M);
	dumsum		= matrix1D(N);
	W 			= matrix2D(N, N);
	wp 			= matrix1D(N);

	for(i = 0; i < N; ++i)
		for(j = 0; j < N; ++j)
			W[i][j] = (double)rand()/((double)RAND_MAX);

	/*FastICA algorithm*/
	for(p = 0; p < N; ++p)
	{
		for(i = 0; i < N; ++i)
			wp[i] = (double)rand()/(double)RAND_MAX;

		VectorNormalization(wp, N);

		/*FastICA iterations*/
		for(it = 0; it < iterations; ++it)
		{
			G = VecMatMult(wp, N, Z, M);

			for(i = 0; i < M; ++i)
				Gder[i] = 1 - tanh(G[i])*tanh(G[i]);

			for(i = 0; i < M; ++i)
				G[i] = tanh(G[i]);

			/*wp = 1/M*Z*G' - 1/M*Gder*ones(M,1)*wp */
			ZGt = MatVecMult(Z, N, M, G);

			for(i = 0; i < N; ++i)
				ZGt[i] /= M;

			GderOnes = 0.0;
			for(i = 0; i < M; ++i)
				GderOnes += Gder[i]/M;

			for(i = 0; i < N; ++i)
				wp[i] = ZGt[i] - GderOnes*wp[i];

			/*Gram-Schmidt decorrelation*/
			for(i = 0; i < N; ++i)
				dumsum[i] = 0.0;

			for(i = 0; i < N; ++i)
			{
				for(j = 0; j < p ; ++j)
				{
					f = 0.0;
					for(k = 0; k < N; ++k)
						f += wp[k]*W[k][j];

					dumsum[i] += f*W[i][j];
				}
			}

			for(i = 0;  i < N; ++i)
				wp[i] -= dumsum[i];

			VectorNormalization(wp, N);
		}

		/*Storing estimated rows of the inverse of the mixing matrix as columns in W*/
		for(i = 0; i < N; ++i)
			W[i][p] = wp[i];
	}

	/*Normalizing estimated inverse of mixing matrix A*/
	for(i = 0; i < N; ++i)
		for(j = 0; j < N; ++j)
			W[i][j] /= sqrt(2.0);

	/*Cleaning up memory allocated, not sure if needed cause we're leaving scope soon*/
	free(G);
	free(Gder);
	free(dumsum),
	free(wp);

	return W;
}
