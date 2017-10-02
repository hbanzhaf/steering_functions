// *** File: utilities.cc
// *** Author(s): Th. Fraichard
// *** Last modified on 04 Aug 2008

// Description: 

#include <utilities.h>

// definition globale d'un epsilon proche de 0 pour les calculs approches
double Epsilon = 0.001;

// gestion d'un epsilon proche de 0 pour les calculs approches
void set_epsilon (double epsilon) { Epsilon = epsilon; }
double get_epsilon () { return Epsilon; }
              
// variables globales de stockage de l'heure (chronometre) 
// timeval (time value) est une structure a deux champs: tv_sec (secondes)
// et tv_usec (microsecondes)
struct timeval time_begin, time_end;
struct timezone time_zone;

// calcul du temps ecoule (en microsecondes) entre deux "timevalues"
double elapsed_time (timeval begin, timeval end) {
  return (end.tv_sec - begin.tv_sec)*1000000 + end.tv_usec - begin.tv_usec; }

// demarrage du chronometre
void start_chronometer () {
  gettimeofday(&time_begin, &time_zone);
}

// mesure du temps ecoule (en microsecondes)
double look_at_chronometer () {
  gettimeofday(&time_end, &time_zone);
  return elapsed_time (time_begin, time_end);
}

// distance cartesienne entre deux points
double point_distance (double x1, double y1, double x2, double y2) { 
  return sqrt (pow (x2 - x1, 2) + pow (y2 - y1, 2)); }

// conversions repere monde <-> repere display			  
double world_to_display_x (double x, double unit_size, double scroll_x) {
  return (x * unit_size + scroll_x); }

double world_to_display_y (double y, double unit_size, double scroll_y) {
  return (-(y * unit_size + scroll_y)); }

double display_to_world_x (double x, double unit_size, double scroll_x) {
  return ((x - scroll_x) / unit_size); }

double display_to_world_y (double y, double unit_size, double scroll_y) {
  return (- (scroll_y + y) / unit_size); }

// conversion angle quelconque -> angle dans l'intervalle [0, 2pi[ radians
double twopify (double alpha) { 
  while ( alpha >= Two_Pi ) { alpha = alpha - (2 * Pi); }
  while ( alpha < 0 ) { alpha = alpha + (2 * Pi); }
  return alpha; }

// conversion radians <-> degres
double degrees (double alpha) { return alpha * 180 / Pi; }
double radians (double alpha) { return alpha * Pi / 180; }

// integrales de fresnel
double fresnel (double x, bool fresnelc) {
#include <fresnel.data>
  if (fabs (x) > Fresnel_Length)
    { cerr << "Fresnel integral out of range" << endl; return 0; }
  int sign = sgn (x);
  double diter = fabs (x) / Fresnel_Length * Fresnel_Samples;
  int iter = (int) diter;
  double rest = diter - iter;
  double inf, sup;
  if ( fresnelc ) { inf = FresnelC[iter]; sup = FresnelC[iter + 1]; }
  else { inf = FresnelS[iter]; sup = FresnelS[iter + 1]; }
  // interpolation
  return (sign * ((1 - rest) * inf + rest * sup)); }

double fresnelc (double s) { return fresnel(s, true); }

double fresnels (double s) { return fresnel(s, false); }

// calcul de l'extremite d'un arc de clothoide
void end_of_clothoid (double x_i, double y_i, 
		      double theta_i, double kappa_i,
		      double sigma, bool forward, double length,
		      double *x_f, double *y_f, 
		      double *theta_f, double *kappa_f) {
  double x, y, theta, kappa, direction;
  bool left;
  if ( forward ) { direction = 1; } else { direction = -1; }
  if ( sigma >= 0 ) { left = true; } else { left = false; }
  // configuration de depart supposee a l'origine (x, y et theta nuls)
  if ( fabs (sigma) < Epsilon ) // sigma nul
    if ( fabs (kappa_i) < Epsilon ) { // et kappa_i nul -> ligne droite
      x = direction * length;
      y = 0;
      theta = 0; 
      kappa = 0;
    }
    else {			// -> arc de cercle
      x = direction * 1 / kappa_i * sin (kappa_i * length);
      y = 1 / kappa_i * (1 - cos (kappa_i * length));
      theta = direction * kappa_i * length;
      kappa = kappa_i;
    }
  else {				// sigma non nul -> clothoide
    // kappa (s) = kappa (0) + sigma * s
    kappa = kappa_i + sigma * length;
    // theta (s) = 1/2*sigma*s^2 + kappa (0)*s
    theta = 0.5 * sigma * pow (length, 2) + kappa_i * length;
    // x (s), y (s): thank you maple...
    // x (s) = int(cos(1/2*sigma*u^2 + kappa*u), u=0..s);
    // y (s) = int(sin(1/2*sigma*u^2 + kappa*u), u=0..s);
    int ssigma = sgn (sigma);	// signe de sigma
    double usigma = fabs (sigma); // unsigned sigma
    int skappa = sgn (kappa_i);	// signe de kappa
    double ukappa = fabs (kappa_i); // unsigned kappa
    double k1 = 0.5 * pow (ukappa, 2) / usigma;
    double k2 = (usigma * length + ssigma * skappa * ukappa) 
      / sqrt (Pi * usigma);
    double k3 = ukappa / sqrt (Pi * usigma);
    x = sqrt (Pi / usigma) * (cos (k1) * fresnelc (k2)
			      + sin (k1) * fresnels (k2) 
			      - ssigma * skappa * cos (k1) * fresnelc (k3) 
			      - ssigma * skappa * sin (k1) * fresnels (k3));
    y = sqrt (Pi / usigma) * (ssigma * cos (k1) * fresnels (k2) 
			      - ssigma * sin (k1) * fresnelc (k2) 
			      - skappa * cos (k1) * fresnels (k3) 
			      + skappa * sin (k1) * fresnelc (k3));
    x = direction * x; 
    theta = direction * theta;
  }
  // translation et rotation pour tenir compte de la configuration de depart
  global_frame_change (x_i, y_i, theta_i, x, y, x_f, y_f);
  *theta_f = twopify (theta_i + theta);
  *kappa_f = kappa;
}

// calcul de l'extremite d'un arc de cercle
void end_of_circular_arc (double x_i, double y_i, double theta_i,
			  double kappa, bool forward, double length,
			  double *x_f, double *y_f, double *theta_f) {
  double x, y, theta, direction;
  if ( forward ) { direction = 1; } else { direction = -1; }
  // configuration de depart supposee a l'origine (x, y et theta nuls)
  if ( fabs (kappa) < Epsilon ) {// kappa nul -> ligne droite
    x = direction * length;
    y = 0;
    theta = 0; 
  }
  else {			// -> arc de cercle
    x = direction * 1 / kappa * sin (kappa * length);
    y = 1 / kappa * (1 - cos (kappa * length));
    theta = direction * kappa * length;
  }
  // translation et rotation pour tenir compte de la configuration de depart
  *x_f = x * cos (theta_i) - y * sin (theta_i) + x_i;
  *y_f = x * sin (theta_i) + y * cos (theta_i) + y_i;
  *theta_f = twopify (theta_i + theta);
}

// changement de repere: soit un repere local defini par son origine (x, y)
// et son orientation theta dans un repere global.  soit (local_x, local_y)
// les coordonnees d'un point dans le repere local, la fonction retourne
// (global_x, global_y) les coordonnees de ce point dans le repere global.
void global_frame_change (double x, double y, double theta, 
			  double local_x, double local_y, 
			  double *global_x, double *global_y) {
  *global_x = local_x * cos (theta) - local_y * sin (theta) + x;
  *global_y = local_x * sin (theta) + local_y * cos (theta) + y; }

// changement de repere: soit un repere local defini par son origine (x, y)
// et son orientation theta dans un repere global.  soit (global_x, global_y)
// les coordonnees d'un point dans le repere global, la fonction retourne
// (local_x, local_y) les coordonnees de ce point dans le repere local.
void local_frame_change (double x, double y, double theta,
			 double global_x, double global_y,
			 double *local_x, double *local_y) {
  *local_x = 0;
  *local_y = 0; }

// "plot" d'une fonction f dans un fichier destine a gnuplot.  f est de
// type "double f (double x)", on plot f de min a max avec le pas
// d'echantillonnage step dans le fichier filename.
void plot_function (double (*f)(double), 
		    double min, double max, double step, 
		    char *filename) {
  ofstream file (filename, ios::out);
  double x = min; double y;
  while ( x <= max ) 
    { 
      y = f (x); 
      file << printf ("%f", x) << " " << printf ("%f", y) << endl; 
      x = x + step; }
  file.close ();
}

// calul de la fonction D1 pour un angle alpha compris entre 0 et pi
double D1 (double alpha)
{ return cos (alpha) * fresnelc (sqrt (2 * alpha / Pi)) + 
    sin (alpha) * fresnels (sqrt (2 * alpha / Pi)); }

// test de la nullite d'un pointeur (p == 0 ?)
bool null (void *p) { return p == 0; }

// determination de l'index du plus petit element d'un tableau de double
int array_index_min (double array [], int size) {
  double min = array [0];
  int index_min = 0;
  for (int i = 1; i < size; i++)
    { if ( array [i] < min ) { index_min = i; min = array [i]; } }
  return index_min;
}

// initialisation a une valeur donnee d'un tableau de double
void double_array_init (double array [], int size, double value) {
  for (int i = 0; i < size; i++) { array [i] = value; }
}

// initialisation a NIL d'un tableau de pointeurs
void pointer_array_init (void * array [], int size) {
  for (int i = 0; i < size; i++) { array [i] = 0; }
}

// initialisation du generateur de nombres aleatoires
void random_init () {
  struct timeval time;
  struct timezone time_zone;
  gettimeofday(&time, &time_zone);
  srand (time.tv_usec);
}

// generation d'un reel aleatoire compris entre 0 et r
double random_real (double r) {
  double random = rand ();
  return random * r / RAND_MAX;
}

// generation d'un reel aleatoire compris entre r1 et r2
double random_real (double r1, double r2) {
  double random = rand ();
  return random * (r2 - r1) / RAND_MAX + r1;
}

