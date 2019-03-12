/*
 * ParameterDefinition.h
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#ifndef PARAMETERDEFINITION_H_
#define PARAMETERDEFINITION_H_

#include <time.h>

	int 	iterations, N, C , M, p;
	float 	finalTime, initialTime, K, na , ns;
	double 	**Amix, **W, **WT, *timeVector, **S, **Sest, **Xobs;
	double 	**X, **Z;
	double 	periodSource5, periodSource6, avgsource5, avgsource6;
	double 	time_spent;
	clock_t begin, end;

#endif /* PARAMETERDEFINITION_H_ */
