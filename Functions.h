/*
 * Functions.h
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

void printMatrix(int N, int M, double **S);
void printMatrixS(int N, int M, double **S);


void ParameterInput();
void setupVars();

double funcSource1(double x);
double funcSource2(double x);
double funcSource3(double x);
double funcSource4(double x);
double funcSource5(double x);
double funcSource6(double x);

void Initialize();
void SetUpSources();
double** XobsGen(double** Amix, double** S, int N, int M);
void ExportingData();
void FreeMemory();

#endif /* FUNCTIONS_H_ */
