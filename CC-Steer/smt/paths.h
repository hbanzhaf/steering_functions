// *** File: paths.h
// *** Author(s): Th. Fraichard
// *** Last modified on 9 Jul 2008

// Description:

#ifndef CC_PATHS_H
#define CC_PATHS_H

#include <iostream>
using namespace std;
#include <math.h>

#include <object.h>
#include <configuration.h>
#include <utilities.h>

// ############################################################################
// classe representant les cc-circles
// ############################################################################
class SMT_CC_Circle : public virtual SMT_Object {
public:
  // configuration origine
  SMT_Configuration start;
  // courbure maximum (signee), derivee de la courbure (signee)
  double kappa, sigma;
  // ca tourne a gauche, a droite ? en marche avant, marche arriere ?
  bool left, forward;
  // centre du cc-circle, rayon
  double xc, yc, radius;
  // ecart angulaire (positif) entre l'orientation de start et la tangente
  // au cc-circle en start
  double mu;
  // configuration intermediaire atteinte par le 1er arc de clotoide
  SMT_Configuration first;

  // constructeur(s)
  // kappa et sigma non obligatoirement signes
  SMT_CC_Circle (SMT_Configuration _start,
		 bool _left, bool _forward,
		 double _kappa, double _sigma);
  // construction d'un cc-circle "libre", i.e. sans start ni first
  // kappa et sigma non obligatoirement signes
  SMT_CC_Circle (double _xc, double _yc, double _radius,
		 bool _left, bool _forward,
		 double _kappa, double _sigma, double _mu);
  // affichage alphanumerique des attributs
  void print (bool eol);
};


// ############################################################################
// generic class for the paths
// ############################################################################
class SMT_Path : public virtual SMT_Object {
public:
  // path length
  double length;
  // computing time (in microseconds) for the path
  double time;

  // constructor
  SMT_Path (double _length, double _time);
};


// ############################################################################
// cc_dubins paths 
// ############################################################################

// path types: EMPTY, S, L, R, etc. were added to simplify certain
// post-processing operations (display, discretization)
enum cc_dubins_path_type {
  LSL, LSR, RSL, RSR, LR1L, LR2L, RL1R, RL2R,
  S, L, R, LS, SL, SR, RS, EMPTY};
const int nb_cc_dubins_paths = 16;

// cc_dubins path class
class SMT_CC_Dubins_Path : public SMT_Path {
public:
  // start and end configurations
  SMT_Configuration start, end;
  // max. curvature (positive), max. curvature derivative (positive)
  double kappa, sigma;
  // path type
  cc_dubins_path_type type;
  // intermediate configurations
  SMT_Configuration *qi1, *qi2;
  // cc-circles possibly associated with the start, end and intermediate
  // configurations
  SMT_CC_Circle *cstart, *cend, *ci1, *ci2;

  // constructor
  SMT_CC_Dubins_Path (SMT_Configuration _start, SMT_Configuration _end,
		      cc_dubins_path_type _type, double kappa, double sigma,
		      SMT_Configuration *_qi1, SMT_Configuration *_qi2,
		      SMT_CC_Circle *_cstart, SMT_CC_Circle *_cend,
		      SMT_CC_Circle *_ci1,
		      double _length, double _time);
  // destructor
  ~SMT_CC_Dubins_Path ();
  // alphanumeric display.  if EOL then end of line
  void print (bool eol);
};


// ############################################################################
// cc_reeds_shepp paths
// ############################################################################

// path types: Segment, T, etc. were added to simplify certain
// post-processing operations (display, discretization)
enum cc_rs_path_type {
  Empty, Segment, T,
  TT, TcT,
  TST, TcST, TScT, TcScT,
  TTT, TcTT, TTcT, TcTcT,
  TTST, TcTST, TTcST, TTScT, TcTcST, TcTScT, TTcScT, TcTcScT,
  TTSTT, TcTSTT, TTcSTT, TTScTT, TTSTcT,
  TcTcSTT, TTcScTT, TTScTcT, TcTScTT, TcTSTcT, TTcSTcT,
  TcTcScTT, TcTcSTcT, TcTScTcT, TTcScTcT, TcTcScTcT
};
const int nb_cc_rs_paths =  37;

// cc_reeds_shepp path class
class SMT_CC_RS_Path : public SMT_Path {
public:
  // start and end configurations
  SMT_Configuration start, end;
  // max. curvature (positive), max. curvature derivative (positive)
  double kappa, sigma;
  // path type
  cc_rs_path_type type;
  // intermediate configurations
  SMT_Configuration *qi1, *qi2, *qi3, *qi4;
  // cc-circles possibly associated with the start, end and intermediate
  // configurations
  SMT_CC_Circle *cstart, *cend, *ci1, *ci2, *ci3;

  // constructor
  SMT_CC_RS_Path (SMT_Configuration _start, SMT_Configuration _end,
		  cc_rs_path_type _type, double kappa, double sigma,
		  SMT_Configuration *_qi1, SMT_Configuration *_qi2,
		  SMT_Configuration *_qi3, SMT_Configuration *_qi4,
		  SMT_CC_Circle *_cstart, SMT_CC_Circle *_cend,
		  SMT_CC_Circle *_ci1, SMT_CC_Circle *_ci2,
		  double _length, double _time);
  // destructor
  ~SMT_CC_RS_Path ();
  // alphanumeric display.  if EOL then end of line
  void print (bool eol);
};


// ############################################################################
// famille des chemins de type cc_tp
// ############################################################################
// classe pour les chemins de type sc
class SMT_SC_Path : public SMT_Path {
public:
  // configurations de depart et d'arrivee
  SMT_Configuration start;
  // sharpness du chemin (signee)
  double sharpness;
  // marche avant, marche arriere ?
  bool forward;

  // constructeur
  SMT_SC_Path (SMT_Configuration _start,
	       double _sharpness, bool _forward, double _length);
  // affichage alphanumerique des attributs
  void print (bool eol);
};

// classe pour les chemins de type cc_tp
class SMT_CC_TP_Path : public SMT_Path {
public:
  // configurations de depart et d'arrivee
  SMT_Configuration start, end;
  // courbure max (positive), derivee de la courbure max (positive)
  double kappa, sigma;
  // configurations intermediaires
  SMT_Configuration *qi1, *qi2;
  // composants de type sc 
  SMT_SC_Path *sci1, *sci2;
  // marche avant, marche arriere entre qi1 et qi2 ?
  bool forward;

  // constructeur
  SMT_CC_TP_Path (SMT_Configuration _start, SMT_Configuration _end,
		  double _kappa, double _sigma,
		  SMT_Configuration *_qi1, SMT_Configuration *_qi2,
		  SMT_SC_Path *_sci1, SMT_SC_Path *_sci2,
		  bool forward, double _length, double _time);
  // destructeur
  ~SMT_CC_TP_Path ();
  // affichage alphanumerique des attributs
  void print (bool eol);
};


// ############################################################################
// fonctions non membres
// ############################################################################

// existence d'une mu-tangente externe entre deux cc-circles
bool external_mu_tangent_exists (SMT_CC_Circle c1, SMT_CC_Circle c2);

// calcul de la mu-tangente externe du cc-circle c1 au cc-circle c2
// (conditions d'existence de la mu-tangente externe supposees etre vraies)
void external_mu_tangent (SMT_CC_Circle c1, SMT_CC_Circle c2,
			  SMT_Configuration **q1, SMT_Configuration **q2);

// existence d'une mu-tangente interne entre deux cc-circles
bool internal_mu_tangent_exists (SMT_CC_Circle c1, SMT_CC_Circle c2);

// calcul de la mu-tangente interne du cc-circle c1 au cc-circle c2
// (conditions d'existence de la mu-tangente interne supposees etre vraies)
void internal_mu_tangent (SMT_CC_Circle c1, SMT_CC_Circle c2,
			  SMT_Configuration **q1, SMT_Configuration **q2);

// existence d'un cc-circle tangent entre deux cc-circles
bool tangent_circle_exists (SMT_CC_Circle c1, SMT_CC_Circle c2);

// calcul des cc-circles tangents entre deux cc-circles
// (conditions d'existence des cc-circles tangents supposees etre vraies)
void tangent_circle (SMT_CC_Circle c1, SMT_CC_Circle c2,
		     SMT_Configuration **q1, SMT_Configuration **q2,
		     SMT_Configuration **q3, SMT_Configuration **q4);

// calcul de la longueur du cc-turn defini par le cc-circle c et la
// configuration d'arrivee q (supposee appartenir au cc-circle et avoir une
// orientation correcte).  le booleen cusps indique si l'on se trouve dans
// le cas avec ou sans manoeuvres (i.e. reeds-shepp vs. dubins)
double cc_turn_length (SMT_CC_Circle *c, SMT_Configuration *q, bool cusps);

// calcul de la longueur du cc-turn de reeds et shepp defini par le
// cc-circle c et la configuration d'arrivee q (supposee appartenir au
// cc-circle et avoir une orientation correcte)
double rs_turn_length (SMT_CC_Circle *c, SMT_Configuration *q);

// calcul de la longueur du cc-turn de dubins defini par le cc-circle c et
// la configuration d'arrivee q (supposee appartenir au cc-circle et avoir
// une orientation correcte)
double dubins_turn_length (SMT_CC_Circle *c, SMT_Configuration *q);


// calcul pour le cc-circle c de la configuration atteinte par un cc-turn
// d'une certaine deflection comprise entre 0 et 2*pi
void deflection_to_configuration (SMT_CC_Circle c, double deflection,
				  SMT_Configuration *q);

// appartenance d'une configuration a un cc-circle
bool configuration_on_cc_circle (SMT_CC_Circle c, 
				 SMT_Configuration q);

#endif
