// *** File: test-nw.cc
// *** Author(s): Th. Fraichard
// *** Last modified on 9 Jul 2008

// Description:

#include <iostream>
#include <fstream>
using namespace std;

#include <utilities.h>
#include <configuration.h>
#include <paths.h>
#include <cc-dubins.h>
#include <cc-rs.h>

// maximum curvature and derivative of curvature
double KappaMax = 1.0; //0.2;
double SigmaMax = 1.0; //0.05;

// start configuration
SMT_Configuration start (0, 0, 0, 0);
// end configuration
SMT_Configuration goal (1, 1, 0.785, 0);


// ############################################################################
int main(int argc, char **argv)
{
  char *message;

  // set goal configuration
  cout << "***** CC-Steer *****\n";
  printf ("Start Config. = (%.2f, %.2f, %.2f, %.2f)\n", 
	  start.x, start.y, start.theta, start.kappa);
  cout << "Enter Goal Config.\n";
  cout << "x = ";
  cin >> goal.x;
  cout << "y = ";
  cin >> goal.y;
  cout << "theta = ";
  cin >> goal.theta;
  cout << "kappa = ";
  cin >> goal.kappa;
  printf ("Goal Config. = (%.2f, %.2f, %.2f, %.2f)\n", 
	  goal.x, goal.y, goal.theta, goal.kappa);

  // compute cc-dubins path
  SMT_CC_Dubins_Path *d_path;
  
  d_path = cc_dubins (start, goal, KappaMax, SigmaMax);
  d_path->print (true);		// alphanumeric display

  // compute cc-reeds-shepp path
  SMT_CC_RS_Path *rs_path;

  rs_path = cc_reeds_shepp (start, goal, KappaMax, SigmaMax);
  rs_path->print (true);	// alphanumeric display
}
