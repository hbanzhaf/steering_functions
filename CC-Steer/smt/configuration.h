// *** File: configuration.h
// *** Author(s): Th. Fraichard
// *** Last modified on 9 Jul 2008

// Description: definition of the configuration class with four parametres:
// position (x,y), orientation and curvature.

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <iostream>
using namespace std;

#include <object.h>
#include <utilities.h>

// ############################################################################
// class SMT_Configuration
// ############################################################################

class SMT_Configuration : public virtual SMT_Object {
public:
  // position, orientation (in radians, between 0 and 2pi) and curvature
  double x, y, theta, kappa;
  // constructor
  SMT_Configuration (double _x = 0.0, double _y = 0.0,
		     double _theta = 0.0, double _kappa = 0.0);
  // are two configurations equal?
  bool operator== (SMT_Configuration q);
  // alphanumeric display.  if EOL then end of line
  void print (bool eol);
};


// ############################################################################
// non members functions
// ############################################################################

// cartesian distance between two configurations
double configuration_distance (SMT_Configuration q1, SMT_Configuration q2);

// are two configurations almost aligned?
bool configuration_aligned (SMT_Configuration q1,
			    SMT_Configuration q2);
// are two configurations almost equal?
bool configuration_equal (SMT_Configuration q1, SMT_Configuration q2);

#endif
