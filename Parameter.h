/*
 * Parameter.h
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#ifndef PARAMETER_H_
#define PARAMETER_H_

#include <time.h>

	extern int 		iterations, N, C , M, p;
	extern float 	finalTime, initialTime, K, na , ns;
	extern double 	**Amix, **W, **WT, *timeVector, **S, **Sest, **Xobs;
	extern double 	**X, **Z;
	extern double 	periodSource5, periodSource6, avgsource5, avgsource6;
	extern double 	time_spent;
	extern clock_t 	begin, end;

#endif /* PARAMETER_H_ */
