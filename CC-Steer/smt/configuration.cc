// *** File: configuration.cc
// *** Author(s): Th. Fraichard
// *** Last modified on 9 Jul 2008

// Description: definition de la classe configuration a quatre parametres:
// position (x,y), orientation et courbure.

#include <configuration.h>

// ############################################################################
// class SMT_Configuration
// ############################################################################

// constructeur(s)
SMT_Configuration::SMT_Configuration (double _x, double _y,
				      double _theta, double _kappa) {
  x = _x; y = _y; theta = twopify (_theta); kappa = _kappa; 
}

// egalite entre deux configurations
bool SMT_Configuration::operator== (SMT_Configuration q) {
  return (x == q.x) && (y == q.y) && (theta == q.theta) && (kappa == q.kappa);
}

// affichage des attributs
void SMT_Configuration::print (bool eol) { 
  cout << "(" << x << ", " << y << ", " 
       << degrees (theta) << ", " << kappa << ")"; 
  if ( eol ) { cout << endl; }
}


// ############################################################################
// fonctions non membres
// ############################################################################

// distance cartesienne entre deux configurations
double configuration_distance (SMT_Configuration q1, 
			       SMT_Configuration q2) {
  return point_distance (q1.x, q1.y, q2.x, q2.y);
}

// deux configurations sont elles alignees ?
bool configuration_aligned (SMT_Configuration q1, 
			    SMT_Configuration q2) {
  // orientations identiques ?
  if ( fabs (q2.theta - q1.theta) > get_epsilon () ) { return false; }
  // angle de la droite orientee joignant q1 a q2
  double angle = twopify (atan2 (q2.y - q1.y, q2.x - q1.x));
  return fabs (angle - q1.theta) <= get_epsilon ();
}

// deux configurations sont-elles egales ?
bool configuration_equal (SMT_Configuration q1, 
			  SMT_Configuration q2) {
  if ( configuration_distance (q1, q2) < get_epsilon () )
    if ( fabs (q2.theta - q1.theta) < get_epsilon () ) { return true; }
  return false;
}
