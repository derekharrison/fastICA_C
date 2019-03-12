/*
 * ParameterInput.c
 *
 *  Created on: 14 apr. 2016
 *      Author: dharrison
 */

#include "../inc/ParameterDefinition.h"


void ParameterInput()
{
    /*
     * The user parameters are entered here. This is the only section of the code
     * Which should be changed by the user. Note that 6 sources are available. If more sources
     * are required or desired these need to be entered in the corresponding function SetUpSources()
     * in Functions.c. The extra sources need to be added in Sources.c
     * If less are required or desired, the excess sources need to be commented or removed
     * from the SetUpSources() function.
     */

    N             = 6;                    //The number of sources. (It is preferable not to change this value)
    C             = N;                    //The number of observations (sample size). This value should be equal to N and should not be changed!
    M             = 10000;                //The number of observation samples

    K             = 0.1;                  //The slope of the zig-zag source
    na            = 8;                    //The amount of zig-zag source periods (amount of peaks)
    ns            = 5;                    //The amount of alternating step-function periods

    finalTime     = 40*3.1415926535;      //Final time (s)
    initialTime   = 0.0;                  //Initial time (s)

    iterations   = 100;                    //Number of FastICA iterations
}
