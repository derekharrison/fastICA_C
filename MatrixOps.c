/*
 * MatrixOps.c
 *
 *  Created on: 15 apr. 2016
 *      Author: dharrison
 */

#include "math.h"
#include "Memory.h"
#include <stdlib.h>


double** MatMult(double** A, int rows1, int columns1, double** B, int rows2, int columns2)
{
	/*
	 * This function performs matrix multiplication of two matrices A and B (A*B)
	 * with at least two rows and columns.
	 * The result Sp = A*B is returned by the function.
	 */

	int i, j, k;
	double **Sp;

	Sp = matrix2D(rows1, columns2);

	for(i = 0; i < rows1; ++i)
		for(j = 0; j < columns2; ++j)
		{
			for(k = 0; k < columns1; ++k)
				Sp[i][j] += A[i][k]*B[k][j];
		}

	return Sp;
}


double* VecMatMult(double* V, int SizeVec, double** B, int columns)
{
	/*
	 * This function performs matrix multiplication of vector V and matrix B (V*B)
	 * with the matrix containing at least two rows and columns
	 * The result Sp = V*B is returned by the function.
	 */

	int i, k;
	double *Sp;

	Sp = matrix1D(columns);

	for(i = 0; i < columns; ++i)
		for(k = 0; k < SizeVec; ++k)
			Sp[i] += V[k]*B[k][i];

	return Sp;
}


double* MatVecMult(double** B, int rows, int columns, double* V)
{
	/*
	 * This function performs matrix multiplication of matrix B and vector V (B*V)
	 * with the matrix containing at least two rows and columns.
	 * The result Sp = B*V is returned by the function.
	 */

	int i, k;
	double *Sp;

	Sp = matrix1D(rows);

	for(i = 0; i < rows; ++i)
		for(k = 0; k < columns; ++k)
			Sp[i] += B[i][k]*V[k];

	return Sp;
}


double** MatTranspose(double** A, int rows, int columns)
{
	/*
	 * This function computes the transpose of matrix A
	 * The tranpose Sp = A' (matlab notation) is returned by the function.
	 */

	int i, j;
	double **Sp;
	Sp = matrix2D(columns, rows);

	for(i = 0; i < columns; ++i)
		for(j = 0; j < rows; ++j)
			Sp[i][j] = A[j][i];

	return Sp;
}


void VectorNormalization(double *wp, int sizeVec)
{
	/*
	 * This function normalizes a vector wp and returns
	 * the normalized vector wp.
	 */

	int i;
	double sqrtwpwp = 0.0;

	for(i = 0; i < sizeVec; ++i)
		sqrtwpwp += wp[i]*wp[i];

	for(i = 0; i < sizeVec; ++i)
		wp[i] = wp[i]/sqrt(sqrtwpwp);
}


void EigenDecomposition(double** ExxT, int N, double** EigVectors, double* EigValues, int iterations)
{
	/*
	 * This function computes the eigenvalues and eigenvectors
	 * of a real, symmetric, N x N matrix.
	 * The eigenvalues are stored in the function argument EigValues and the eigenvectors are
	 * stored in the function argument EigVectors.
	 */

	int i, j, k, p, it;
	double **EigVecs, **Q, **EigVals, *wp, *dumsum;
	double **R, **QT;
	double f;

	EigVecs		= matrix2D(N, N);
	Q 			= matrix2D(N, N);
	EigVals 	= matrix2D(N, N);
	wp 			= matrix1D(N);
	dumsum		= matrix1D(N);
	R			= matrix2D(N, N);
	QT			= matrix2D(N, N);

	/*Initializing Ait, matrix containing eigenvalues as the diagonal*/
	for(i = 0; i < N; ++i)
		for(j = 0; j < N; ++j)
			EigVals[i][j] = ExxT[i][j];

	/*Initializing Q and E for computation of eigenvectors*/
	for(i = 0; i < N; ++i)
	{
		Q[i][i] = 1.0;
		EigVecs[i][i] = 1.0;
	}

	/*Eigen decomposition iterations*/
	for(it = 0; it < iterations; ++it)
	{
		/*Gram-Schmidt decorrelation*/
		for(p = 0; p < N; ++p)
		{
			for(i = 0; i < N; ++i)
				wp[i] = EigVals[i][p];

			VectorNormalization(wp, N);

			for(i = 0; i < N; ++i)
				dumsum[i] = 0.0;

			for(i = 0; i < N; ++i)
			{
				for(j = 0; j < p ; ++j)
				{
					f = 0.0;
					for(k = 0; k < N; ++k)
						f += wp[k]*Q[k][j];

					dumsum[i] += f*Q[i][j];
				}
			}

			for(i = 0;  i < N; ++i)
				wp[i] -= dumsum[i];

			VectorNormalization(wp, N);

			/*Storing estimated rows of the inverse of the mixing matrix as columns in W*/
			for(i = 0; i < N; ++i)
				Q[i][p] = wp[i];
		}

		QT = MatTranspose(Q, N, N);

		R = MatMult(QT, N, N, EigVals, N, N);

		EigVals = MatMult(R, N, N, Q, N, N);

		EigVecs = MatMult(EigVecs, N, N, Q, N, N);
	}

	EigVecs = MatMult(EigVecs, N, N, Q, N, N);

	for(i = 0; i < N; ++i)
		EigValues[i] = EigVals[i][i];

	for(i = 0; i < N; ++i)
		for(j = 0; j < N; ++j)
			EigVectors[i][j] = EigVecs[i][j];

	free(Q);
	free(EigVals);
	free(wp);
	free(dumsum);
	free(R);
	free(QT);
	free(EigVecs);
}
