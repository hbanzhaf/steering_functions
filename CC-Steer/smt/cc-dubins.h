// *** File: cc-dubins.h
// *** Author(s): Th. Fraichard
// *** Last modified on 06 Aug 2008

// Description: 

#ifndef CC_DUBINS_H
#define CC_DUBINS_H

#include <iostream>
using namespace std;

#if __GNUC__ == 4
#include <values.h>	      /* to have MAXFLOAT, FLT_MAX with gcc 4.x */
#endif
#if __GNUC__ == 3
#include <float.h>	      /* to have MAXFLOAT, FLT_MAX with gcc 3.x */
#endif

#include <utilities.h>
#include <configuration.h>
#include <paths.h>

SMT_CC_Dubins_Path *cc_dubins (SMT_Configuration start, 
			       SMT_Configuration end,
			       double kappa, double sigma);

#endif
