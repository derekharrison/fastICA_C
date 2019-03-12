/*
 * ExportingData.c
 *
 *  Created on: 15 apr. 2016
 *      Author: dharrison
 */

#include "Parameter.h"
#include <stdio.h>

void ExportingData()
{
	/*
	 * This function exports the data obtained from the FastICA algorithm
	 * to a textfile 'SourcesEstimation.txt'.
	 * The data in the textfile is visualized in Matlab using the Matlab script
	 * ReadingData.m
	 */

	int i, j;

	/*Exporting data*/
	  FILE *pfile;

	  pfile = fopen("SourcesEstimation.txt","w");
	  if (pfile != NULL)
	  {
		 fprintf(pfile,"%d\n%d\n",N, M);
  		 for(i = 0; i < M; ++i)
  		 {
  			 for(j = 0; j < N; ++j)
  				 fprintf(pfile,"%.2f\t", Sest[j][i]);
	  		 fprintf(pfile,"\n");
  		 }

  		 for(i = 0; i < M; ++i)
  			fprintf(pfile,"%.2f\n", timeVector[i]);

		 fclose(pfile);
	  }
	  else
	  {
		 printf("Could not open file");
	  }
	  printf("\n\n");
}
