/*
 * Sources.c
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#include <math.h>
#include "Parameter.h"


double funcSource1(double x)
{
	return sin(1.1*x);
}


double funcSource2(double x)
{
	return cos(0.25*x);
}


double funcSource3(double x)
{
	return sin(0.1*x);
}


double funcSource4(double x)
{
	return cos(0.7*x);
}


double funcSource5(double x)
{
	return K*x - floor(x/periodSource5)*K*periodSource5;
}


double funcSource6(double x)
{
	if((int)floor(x/periodSource6) % 2 == 0)
		return 1;
	else
		return -1;
}
