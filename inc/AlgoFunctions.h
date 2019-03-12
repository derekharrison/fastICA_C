/*
 * AlgoFunctions.h
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#ifndef ALGOFUNCTIONS_H_
#define ALGOFUNCTIONS_H_

double** PreprocessingCentering(double** Xobs, int rows, int columns);
double** PreprocessingWhitening(double** X, int N, int M);
double** FastICA(double ** Z, int N, int M, int iterations);

#endif /* ALGOFUNCTIONS_H_ */
