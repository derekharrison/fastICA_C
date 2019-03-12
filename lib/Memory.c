/*
 * Memory.c
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#include <stdlib.h>


double *matrix1D(int np)
{
    /*
     * This function allocates memory for a vector of doubles.
     * A pointer a to the memory allocated is returned by the
     * function.
     */

    double *a;

    a = (double *) calloc(np, sizeof(double));

    return a;
}


double **matrix2D(int nm, int np)
{
    /*
     * This function allocates memory for a matrix of doubles.
     * A pointer m to the memory allocated is returned by the
     * function.
     */

   int i;
   double **m;

   m = (double **) calloc ( nm, sizeof( double *));
   for ( i = 0; i < nm; i++)
      m[i] = (double *) calloc ( np, sizeof( double));

   return m;
}
