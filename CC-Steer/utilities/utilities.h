// *** File: utilities.h
// *** Author(s): Th. Fraichard
// *** Last modified on 04 Aug 2008

// Description: 

#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <fstream>
using namespace std;

#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

/// set epsilon's value (for numerical approximation).
void set_epsilon (double epsilon);
/// get epsilon's value (for numerical approximation).
double get_epsilon ();
              
/// constant approximation of Pi.
#define	Pi      3.1415926535897932384
/// constant approximation of Pi/2.
#define	Half_Pi 1.5707963267948966192
/// constant approximation of 2*Pi.
#define	Two_Pi  6.2831853071795864770
/// constant approximation of sqrt(Pi).
#define Sqrt_Pi 1.7724538509055160273

/// sign of a number: -1 if negatve, 1 if positive or zero.
#define sgn(x) (((x) < 0) ? -1 : 1)

/// cartesian distance between two points.
double point_distance (double x1, double y1, double x2, double y2);

/// start chronometer.
void start_chronometer ();

/// get elasped time (in microseconds) since the chronometer was last started.
double look_at_chronometer ();

/// \brief conversions repere "monde" <-> repere "display".
///
/// conversions repere "monde" <-> repere "display".  une unite est
/// representee par unit_size pixels.  le repere "monde" se trouve a la
/// position (scroll_x, scroll_y) par rapport au repere "display".
double world_to_display_x (double x, double unit_size, double scroll_x);
double world_to_display_y (double y, double unit_size, double scroll_y);
double display_to_world_x (double x, double unit_size, double scroll_x); 
double display_to_world_y (double y, double unit_size, double scroll_y);

/// conversion angle quelconque -> angle dans l'intervalle [0, 2pi[ radians
double twopify (double alpha);

/// conversion radians <-> degres
double degrees (double alpha);
double radians (double alpha);

/// integrales de fresnel
double fresnelc (double s);
double fresnels (double s);

/// calcul de l'extremite d'un arc de clotoide
/// x_i, y_i, theta_i, kappa_i: configuration initiale
/// sigma: "sharpness" de la clotoide (signee)
/// forward: marche avant, marche arriere ?
/// length: longueur de l'arc de clotoide (positive)
/// resultat retourne dans x_f, y_f, theta_f, kappa_f
void end_of_clothoid (double x_i, double y_i, 
		      double theta_i, double kappa_i,
		      double sigma, bool forward, double length,
		      double *x_f, double *y_f, 
		      double *theta_f, double *kappa_f);

/// calcul de l'extremite d'un arc de cercle
/// x_i, y_i, theta_i: configuration initiale
/// kappa: courbure du cercle (signee)
/// forward: marche avant, marche arriere ?
/// length: longueur de l'arc de cercle (positive)
/// resultat retourne dans x_f, y_f, theta_f
void end_of_circular_arc (double x_i, double y_i, double theta_i,
			  double kappa, bool forward, double length,
			  double *x_f, double *y_f, double *theta_f);

/// changement de repere: soit un repere local defini par son origine (x, y)
/// et son orientation theta dans un repere global.  soit (local_x, local_y)
/// les coordonnees d'un point dans le repere local, la fonction retourne
/// (global_x, global_y) les coordonnees de ce point dans le repere global.
void global_frame_change (double x, double y, double theta, 
			  double local_x, double local_y, 
			  double *global_x, double *global_y);

/// changement de repere: soit un repere local defini par son origine (x, y)
/// et son orientation theta dans un repere global.  soit (global_x, global_y)
/// les coordonnees d'un point dans le repere global, la fonction retourne
/// (local_x, local_y) les coordonnees de ce point dans le repere local.
void local_frame_change (double x, double y, double theta,
			 double global_x, double global_y,
			 double *local_x, double *local_y);

/// "plot" d'une fonction f dans un fichier destine a gnuplot.  f est de
/// type "double f (double x)", on plot f de min a max avec le pas
/// d'echantillonnage step dans le fichier filename.
void plot_function (double (*f)(double),
		    double min, double max, double step,
		    char *filename);

/// calul de la fonction D1 pour un angle alpha compris entre 0 et pi
double D1 (double alpha);

/// test de la nullite d'un pointeur (p == 0 ?)
bool null (void *p);

/// determination de l'index du plus petit element d'un tableau de double
int array_index_min (double array [], int size);

/// initialisation a une valeur donnee d'un tableau de double
void double_array_init (double array [], int size, double value);

/// initialisation a NIL d'un tableau de pointeurs
void pointer_array_init (void * array [], int size);

/// initialisation du generateur de nombres aleatoires
void random_init ();

/// generation d'un reel aleatoire compris entre 0 et r
double random_real (double r);

/// generation d'un reel aleatoire compris entre r1 et r2
double random_real (double r1, double r2);

#endif
