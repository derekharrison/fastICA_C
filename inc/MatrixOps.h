/*
 * MatrixOps.h
 *
 *  Created on: 15 apr. 2016
 *      Author: dharrison
 */

#ifndef MATRIXOPS_H_
#define MATRIXOPS_H_

double** MatTranspose(double** A, int rows, int columns);
void VectorNormalization(double *wp, int sizeVec);
double* VecMatMult(double* V, int SizeVec, double** B, int columns);
double* MatVecMult(double** B, int rows, int columns, double* V);
double** MatMult(double** A, int rows1, int columns1, double** B, int rows2, int columns2);
void EigenDecomposition(double** ExxT, int N, double** EigVectors, double* EigValues, int iterations);

#endif /* MATRIXOPS_H_ */
