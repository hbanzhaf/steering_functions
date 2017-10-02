// *** File: paths.cc
// *** Author(s): Th. Fraichard
// *** Last modified on 3 May 2001

// Description: 

#include <paths.h>

// ###########################################################################
// class SMT_Path
// ###########################################################################
// constructeur(s)
SMT_Path::SMT_Path (double _length, double _time) {
  length = _length; time = _time; }


// ###########################################################################
// class SMT_CC_Dubins_Path
// ###########################################################################
// constructeur(s)
SMT_CC_Dubins_Path::SMT_CC_Dubins_Path (SMT_Configuration _start, 
					SMT_Configuration _end,
					cc_dubins_path_type _type, 
					double _kappa, double _sigma,
					SMT_Configuration *_qi1, 
					SMT_Configuration *_qi2,
					SMT_CC_Circle *_cstart, 
					SMT_CC_Circle *_cend, 
					SMT_CC_Circle *_ci1, 
					double _length, double _time)
  : SMT_Path (_length, _time) {
  start = _start; end = _end; type = _type; kappa = _kappa; sigma = _sigma;
  qi1 = _qi1; qi2 = _qi2; 
  cstart = _cstart; cend = _cend; ci1 = _ci1; ci2 = 0;
}	      

// destructeur
SMT_CC_Dubins_Path::~SMT_CC_Dubins_Path () {
   delete qi1; delete qi2; 
   delete cstart; delete ci1; delete ci2; delete cend;
}

// affichage alphanumerique des attributs
void SMT_CC_Dubins_Path::print (bool eol) {
  cout << "* CC_Dubins_Path: ";
  start.print (false); cout << " -> "; 
  if ( !null (qi1) ) { qi1->print (false); cout << " -> "; }
  if ( !null (qi2) ) { qi2->print (false); cout << " -> "; }
  end.print (false); cout << ", type ";
  switch (type ) { 
  case LSL: cout << "LSL"; break;
  case LSR: cout << "LSR"; break;
  case RSL: cout << "RSL"; break;
  case RSR: cout << "RSR"; break;
  case LR1L: cout << "LR1L"; break;
  case LR2L: cout << "LR2L"; break;
  case RL1R: cout << "RL1R"; break; 
  case RL2R: cout << "RL2R"; break;
  case S: cout << "S"; break;
  case L: cout << "L"; break;
  case R: cout << "R"; break;
  case LS: cout << "LS"; break;
  case SL: cout << "SL"; break;
  case SR: cout << "SR"; break;
  case RS: cout << "RS"; break;
  case EMPTY: cout << "EMPTY"; break;
  default: cout << "?"; break;}
  cout << ", length " << length << ", time " << time; 
  // eventuel saut de ligne final
  if ( eol ) { cout << endl; }
}


// ###########################################################################
// class SMT_CC_Reeds_Shepp_Path
// ###########################################################################
// constructeur(s)
SMT_CC_RS_Path::SMT_CC_RS_Path (SMT_Configuration _start, 
				SMT_Configuration _end,
				cc_rs_path_type _type, 
				double _kappa, double _sigma,
				SMT_Configuration *_qi1, 
				SMT_Configuration *_qi2, 
				SMT_Configuration *_qi3, 
				SMT_Configuration *_qi4,
				SMT_CC_Circle *_cstart, 
				SMT_CC_Circle *_cend, 
				SMT_CC_Circle *_ci1, 
				SMT_CC_Circle *_ci2, 
				double _length, double _time) 
: SMT_Path (_length, _time) {
  start = _start; end = _end; type = _type; kappa = _kappa; sigma = _sigma;
  qi1 = _qi1; qi2 = _qi2; qi3 = _qi3; qi4 = _qi4; 
  cstart = _cstart; cend = _cend; ci1 = _ci1; ci2 = _ci2; ci3 = 0;  
}

// destructeur
SMT_CC_RS_Path::~SMT_CC_RS_Path () {
  delete qi1; delete qi2; delete qi3; delete qi4;
  delete cstart; delete ci1; delete ci2; delete ci3; delete cend;
}

// affichage alphanumerique des attributs
void SMT_CC_RS_Path::print (bool eol) {
  cout << "* CC_RS_Path: type ";
  switch ( type ) { 
  case Empty: cout << "Empty"; break;
  case Segment: cout << "Segment"; break;
  case T: cout << "T"; break;
  case TT: cout << "TT"; break;
  case TcT: cout << "TcT"; break;
  case TST: cout << "TST"; break;
  case TcST: cout << "TcST"; break;
  case TScT: cout << "TScT"; break;
  case TcScT: cout << "TcScT"; break;
  case TTT: cout << "TTT"; break;
  case TcTT: cout << "TcTT"; break;
  case TTcT: cout << "TTcT"; break;
  case TcTcT: cout << "TcTcT"; break;
  case TTST: cout << "TTST"; break;
  case TcTST: cout << "TcTST"; break;
  case TTcST: cout << "TTcST"; break;
  case TTScT: cout << "TTScT"; break;
  case TcTcST: cout << "TcTcST"; break;
  case TcTScT: cout << "TcTScT"; break;
  case TTcScT: cout << "TTcScT"; break;
  case TcTcScT: cout << "TcTcScT"; break;
  case TTSTT: cout << "TTSTT"; break;
  case TcTSTT: cout << "TcTSTT"; break;
  case TTcSTT: cout << "TTcSTT"; break;
  case TTScTT: cout << "TTScTT"; break;
  case TTSTcT: cout << "TTSTcT"; break;
  case TcTcSTT: cout << "TcTcSTT"; break;
  case TTcScTT: cout << "TTcScTT"; break;
  case TTScTcT: cout << "TTScTcT"; break;
  case TcTScTT: cout << "TcTScTT"; break;
  case TcTSTcT: cout << "TcTSTcT"; break;
  case TTcSTcT: cout << "TTcSTcT"; break;
  case TcTcScTT: cout << "TcTcScTT"; break;
  case TcTcSTcT: cout << "TcTcSTcT"; break;
  case TcTScTcT: cout << "TcTScTcT"; break;
  case TTcScTcT: cout << "TTcScTcT"; break;
  case TcTcScTcT: cout << "TcTcScTcT"; break;
  default: cout << "?"; break; }
  //
  cout << ", length " << length << ", time " << time << " - ";
  // configurations de depart et d'arrivee
  start.print (false); cout << " -> "; 
  if ( !null (qi1) ) { qi1->print (false); cout << " -> "; }
  if ( !null (qi2) ) { qi2->print (false); cout << " -> "; }
  if ( !null (qi3) ) { qi3->print (false); cout << " -> "; }
  if ( !null (qi4) ) { qi4->print (false); cout << " -> "; }
  end.print (false); 
  // informations supplementaires en fonction du type
  if ( true ) {
    switch (type ) { 
    case Empty: break;
    case Segment: break;
    case T: 
      cout << endl; 
      cstart->print (true); break;
    case TT: case TcT: case TST: case TcST: case TScT: case TcScT: 
      cout << endl; 
      cstart->print (true); cend->print (false); break;
    case TTT: case TcTT: case TTcT: case TcTcT: 
      cout << endl; 
      cstart->print (true); ci1->print (true); cend->print (false); 
      break;
    case TTST: case TcTST: case TTcST: case TTScT: 
    case TcTcST: case TcTScT: case TTcScT: case TcTcScT: 
      cout << endl; 
      cstart->print (true); ci1->print (true); cend->print (false); 
      break;
    case TTSTT: case TcTSTT: case TTcSTT: case TTScTT: 
    case TTSTcT: case TcTcSTT: case TTcScTT: case TTScTcT: 
    case TcTScTT: case TcTSTcT: case TTcSTcT: case TcTcScTT: 
    case TcTcSTcT: case TcTScTcT: case TTcScTcT: case TcTcScTcT: 
      cout << endl; 
      cstart->print (true); ci1->print (true); 
      ci2->print (true); cend->print (false); 
      break;
    default: break; }
  }
  // eventuel saut de ligne final
  if ( eol ) { cout << endl; }
}

// ###########################################################################
// class SMT_CC_Circle
// ###########################################################################

// constructeur(s)
SMT_CC_Circle::SMT_CC_Circle (SMT_Configuration _start,
			      bool _left, bool _forward,
			      double _kappa, double _sigma) {
  start = _start; start.kappa = 0;
  left = _left; forward = _forward;
  if ( left ) { kappa = fabs (_kappa); sigma = fabs (_sigma); } 
  else { kappa = - fabs (_kappa); sigma = - fabs (_sigma); } 
  // configuration intermediaire atteinte par le 1er arc de clotoide
  double x_i, y_i, theta_i, kappa_i;
  end_of_clothoid (start.x, start.y, start.theta, 0,
		   sigma, forward, fabs (kappa / sigma), 
		   &x_i, &y_i, &theta_i, &kappa_i);
  first = SMT_Configuration (x_i, y_i, theta_i, kappa_i);
  // centre du cc-circle, rayon
  xc = x_i - sin (theta_i) / kappa;
  yc = y_i + cos (theta_i) / kappa;
  radius  = point_distance (xc, yc, start.x, start.y);
  // ecart angulaire entre l'orientation de start et la tangente au
  // cc-circle en start
  // vecteur d'orientation _start.theta ()
  double u_x = cos (start.theta); double u_y = sin (start.theta);
  double u_norm = sqrt (pow (u_x, 2) + pow (u_y, 2));
  // vecteur normal au rayon
  double v_x; double v_y; 
  if ( left )
    { v_x = yc - start.y; v_y = -xc + start.x; }
  else 
    { v_x = - yc + start.y; v_y = xc - start.x; }
  double v_norm = sqrt (pow (v_x, 2) + pow (v_y, 2));
  // produit scalaire -> cos (angle) -> angle
  double dot_product = (u_x * v_x + u_y * v_y) / (u_norm * v_norm);
  mu = acos (dot_product);  
}

// construction d'un cc-circle "libre", i.e. sans start ni first
SMT_CC_Circle::SMT_CC_Circle (double _xc, double _yc, double _radius,
			      bool _left, bool _forward,
			      double _kappa, double _sigma, double _mu) {
  start = SMT_Configuration (0, 0, 0, 0);
  left = _left; forward = _forward;
  if ( left ) { kappa = fabs (_kappa); sigma = fabs (_sigma); } 
  else { kappa = - fabs (_kappa); sigma = - fabs (_sigma); } 
  first = SMT_Configuration (0, 0, 0, 0);
  xc = _xc; yc = _yc; radius = _radius; mu = _mu;
}

// affichage alphanumerique des attributs
void SMT_CC_Circle::print (bool eol) {
  cout << "CC_Circle: ";
  start.print (false); 
  cout << ", kappa: " << kappa << ", sigma: " << sigma;
  if ( left ) { cout << ", left"; } else { cout << ", right"; } 
  if ( forward ) { cout << ", forward, "; } else { cout << ", backward, "; } 
  first.print (false); 
  cout << ", centre: (" << xc << ", " << yc << "), radius " << radius 
       << ", mu: " << degrees (mu); 
  if ( eol ) { cout << endl; }
}

// existence d'une mu-tangente externe entre deux cc-circles
bool external_mu_tangent_exists (SMT_CC_Circle c1, SMT_CC_Circle c2) {
  // les deux cc-circles ont-ils les meme rayons et angles tangents ?
  if ( fabs (c1.radius - c2.radius) > get_epsilon () ) { return false; }
  if ( fabs (c1.mu - c2.mu) > get_epsilon () ) { return false; }
  // les deux cc-circles sont-ils compatibles ?
  if ( c1.left != c2.left ) { return false; }
  if ( c1.forward == c2.forward ) { return false; }
  // les deux centres sont-ils suffisamment eloignes ?
  double distance = point_distance (c1.xc, c1.yc, c2.xc, c2.yc);
  return ( distance >= 2 * c1.radius * sin (c1.mu)); 
}

// calcul de la mu-tangente externe du cc-circle c1 au cc-circle c2
// (conditions d'existence de la mu-tangente externe supposees etre vraies)
void external_mu_tangent (SMT_CC_Circle c1, SMT_CC_Circle c2,
			  SMT_Configuration **q1, SMT_Configuration **q2) {
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2.yc - c1.yc, c2.xc - c1.xc);
  // parametres de position des config. tangentes dans le repere local au
  // centre de c1
  double delta_x = fabs (c1.radius * sin (c1.mu));
  double delta_y = fabs (c1.radius * cos (c1.mu));
  //
  double x, y;
  if ( c1.left && c1.forward )
    {
      global_frame_change (c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q1 = new SMT_Configuration (x, y, theta, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new SMT_Configuration (x, y, theta, 0);
    }
  if ( !c1.left && c1.forward )
    {
      global_frame_change (c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q1 = new SMT_Configuration (x, y, theta, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new SMT_Configuration (x, y, theta, 0);
    }
}

// existence d'une mu-tangente interne entre deux cc-circles
bool internal_mu_tangent_exists (SMT_CC_Circle c1, SMT_CC_Circle c2) {
  // les deux cc-circles sont-ils compatibles ?
  if ( c1.left == c2.left ) { return false; }
  if ( c1.forward == c2.forward ) { return false; }
  // les deux centres sont-ils suffisamment eloignes ?
  double distance = point_distance (c1.xc, c1.yc, c2.xc, c2.yc);
  return ( distance >= 2 * c1.radius); 
}

// calcul de la mu-tangente interne du cc-circle c1 au cc-circle c2
// (conditions d'existence de la mu-tangente interne supposees etre vraies)
void internal_mu_tangent (SMT_CC_Circle c1, SMT_CC_Circle c2,
			  SMT_Configuration **q1, SMT_Configuration **q2) {
  // distance des centres
  double distance = point_distance (c1.xc, c1.yc, c2.xc, c2.yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2.yc - c1.yc, c2.xc - c1.xc);
  // angle entre la droite des centres et la direction de la mu-tangente
  double alpha = fabs (asin (2 * c1.radius * cos (c1.mu) / distance));
  // parametres de position des config. tangentes dans le repere local au
  // centre de c1
  double delta_x = fabs (c1.radius * sin (c1.mu));
  double delta_y = fabs (c1.radius * cos (c1.mu));
  //
  double x, y;
  if ( c1.left && c1.forward )
    {
      theta = theta + alpha;
      global_frame_change (c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q1 = new SMT_Configuration (x, y, theta, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new SMT_Configuration (x, y, theta, 0);
    }
  if ( !c1.left && c1.forward )
    {
      theta = theta - alpha;
      global_frame_change (c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q1 = new SMT_Configuration (x, y, theta, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new SMT_Configuration (x, y, theta, 0);
    }
}

// existence d'un cc-circle tangent entre deux cc-circles
bool tangent_circle_exists (SMT_CC_Circle c1, SMT_CC_Circle c2) {
  // les deux cc-circles sont-ils compatibles ?
  if ( c1.left != c2.left ) { return false; }
  if ( c1.forward == c2.forward ) { return false; }
  // les deux centres sont-ils suffisamment proches ?
  double distance = point_distance (c1.xc, c1.yc, c2.xc, c2.yc);
  return ( distance <= 4 * c1.radius); 
}

// calcul des cc-circles tangents entre deux cc-circles
// (conditions d'existence des cc-circles tangents supposees etre vraies)
void tangent_circle (SMT_CC_Circle c1, SMT_CC_Circle c2,
		     SMT_Configuration **q1, SMT_Configuration **q2,
		     SMT_Configuration **q3, SMT_Configuration **q4) {
  // distance des centres
  double distance = point_distance (c1.xc, c1.yc, c2.xc, c2.yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2.yc - c1.yc, c2.xc - c1.xc);
  // distance du centre du cc-circle a la droite des centres
  double h = sqrt (fabs (pow (2 * c1.radius, 2) - pow (0.5 * distance, 2)));
  // angle entre la droite des centres et le rayon au centre du cc-circle
  double alpha = fabs (atan (2 * h / distance));
  // parametres de position des config. tangentes dans le repere local au
  // centre de c1
  double delta_x = fabs (c1.radius * cos (alpha));
  double delta_y = fabs (c1.radius * sin (alpha));
  //
  double x, y;
  if ( c1.left && c1.forward )
    {
//        cout << "left-forward: " << distance << ", " << degrees (theta)
//  	   << ", " << h << ", " << degrees (alpha) << ", " << delta_x 
//  	   << ", " << delta_y << endl;
      global_frame_change (c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q1 = new SMT_Configuration (x, y, theta + alpha + Half_Pi - c1.mu, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new SMT_Configuration (x, y, theta - alpha + 1.5 * Pi + c1.mu, 0);
      global_frame_change (c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q3 = new SMT_Configuration (x, y, theta - alpha + Half_Pi - c1.mu, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q4 = new SMT_Configuration (x, y, theta + alpha + 1.5 * Pi + c1.mu, 0);
    }
  if ( c1.left && !c1.forward )
    {
      global_frame_change (c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q1 = new SMT_Configuration (x, y, theta - alpha + Half_Pi + c1.mu, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new SMT_Configuration (x, y, theta + alpha + 1.5 * Pi - c1.mu, 0);
      global_frame_change (c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q3 = new SMT_Configuration (x, y, theta + alpha + Half_Pi + c1.mu, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q4 = new SMT_Configuration (x, y, theta - alpha + 1.5 * Pi - c1.mu, 0);
    } 
  if ( !c1.left && c1.forward )
    {
      global_frame_change (c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q1 = new SMT_Configuration (x, y, theta - alpha - Half_Pi + c1.mu, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q2 = new SMT_Configuration (x, y, theta + alpha + Half_Pi - c1.mu, 0);
      global_frame_change (c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q3 = new SMT_Configuration (x, y, theta + alpha - Half_Pi + c1.mu, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q4 = new SMT_Configuration (x, y, theta - alpha + Half_Pi - c1.mu, 0);
    }
  if ( !c1.left && !c1.forward )
    {
      global_frame_change (c1.xc, c1.yc, theta, delta_x, delta_y, &x, &y);
      *q1 = new SMT_Configuration (x, y, theta + alpha - Half_Pi - c1.mu, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, delta_y, &x, &y);
      *q2 = new SMT_Configuration (x, y, theta - alpha + Half_Pi + c1.mu, 0);
      global_frame_change (c1.xc, c1.yc, theta, delta_x, -delta_y, &x, &y);
      *q3 = new SMT_Configuration (x, y, theta - alpha - Half_Pi - c1.mu, 0);
      global_frame_change (c2.xc, c2.yc, theta, -delta_x, -delta_y, &x, &y);
      *q4 = new SMT_Configuration (x, y, theta + alpha - Half_Pi + c1.mu, 0);
    }
}


// ###########################################################################
// class SMT_CC_Turn
// ###########################################################################

// calcul de la longueur du cc-turn defini par le cc-circle c et la
// configuration d'arrivee q (supposee appartenir au cc-circle et avoir une
// orientation correcte).  le booleen cusps indique si l'on se trouve dans
// le cas avec ou sans manoeuvres (i.e. reeds-shepp vs. dubins)
double cc_turn_length (SMT_CC_Circle *c, SMT_Configuration *q, bool cusps) {
  double deflection_min = twopify (pow (c->kappa, 2) / fabs (c->sigma));
  // calcul de la deflection
  double deflection;
  double alpha_c = c->start.theta;
  double alpha_q = q->theta;
  if ( c->left && c->forward ) { deflection = twopify (alpha_q - alpha_c); }
  if ( c->left && !c->forward ) { deflection = twopify (alpha_c - alpha_q); }
  if ( !c->left && c->forward ) { deflection = twopify (alpha_c - alpha_q); }
  if ( !c->left && !c->forward ) { deflection = twopify (alpha_q - alpha_c); }

  // deflection nulle -> cc-turn = segment
  if ( fabs (deflection) < get_epsilon () ) { 
    return fabs (2 * c->radius * sin (c->mu)); }

  // cc-turn irregulier
  if ( deflection < deflection_min ) {
    double d1 = D1 (deflection / 2);
    double d2 = point_distance (c->start.x, c->start.y, 
				q->x, q->y);
    double sharpness = 4 * Pi * pow (d1, 2) / pow (d2, 2);
    return 2 * sqrt (deflection / sharpness);
  }

  // cc-turn regulier (cas sans manoeuvres)
  double length_min = fabs (2 * c->kappa / c->sigma);
  switch ( cusps ) {		// cc-turn regulier
  case true:			// cas avec manoeuvres
    if ( deflection < deflection_min + Pi ) 
      { return length_min + fabs ((deflection - deflection_min) / c->kappa); }
    else
      { return length_min + 
	  fabs ((Two_Pi - deflection + deflection_min) / c->kappa); }    
    break;
  default:			// cas sans manoeuvres
    return length_min + fabs ((deflection - deflection_min) / c->kappa);
    break; }
}

// calcul de la longueur du cc-turn de reeds et shepp defini par le
// cc-circle c et la configuration d'arrivee q (supposee appartenir au
// cc-circle et avoir une orientation correcte)
double rs_turn_length (SMT_CC_Circle *c, SMT_Configuration *q) {
  return cc_turn_length (c, q, true);
}

// calcul de la longueur du cc-turn de dubins defini par le cc-circle c et
// la configuration d'arrivee q (supposee appartenir au cc-circle et avoir
// une orientation correcte)
double dubins_turn_length (SMT_CC_Circle *c, SMT_Configuration *q) {
  return cc_turn_length (c, q, false);
}

// calcul pour le cc-circle c de la configuration atteinte par un cc-turn
// d'une certaine deflection comprise entre 0 et 2*pi
void deflection_to_configuration (SMT_CC_Circle c, double deflection,
				  SMT_Configuration *q) {
  double x, y, theta, x_g, y_g, theta_g;
  // on se place a la configuration origine, on tourne a gauche et en avant 
  x = c.radius * (sin (deflection + c.mu) + sin (c.mu));
  y = c.radius * (cos (-c.mu) - cos (deflection + c.mu));
  theta  = deflection;
  // on tient compte de la veritable nature du virage (gauche / droite,
  // marche avant / arriere)
  if ( c.left && !c.forward ) { x = -x; theta = -theta; }
  if ( !c.left && c.forward ) { y = -y; theta = -theta; }
  if ( !c.left && !c.forward ) { x = -x; y = -y; }
  // on tient compte de la position reelle de le configuration de depart
  global_frame_change (c.start.x, c.start.y, 
		       c.start.theta, 
		       x, y, &x_g, &y_g);
  theta_g = theta + c.start.theta;
  // retour du resultat
  q->x = x_g;
  q->y = y_g;
  q->theta = twopify (theta_g);
  q->kappa = 0;
}

// appartenance d'une configuration a un cc-circle
bool configuration_on_cc_circle (SMT_CC_Circle c, 
				 SMT_Configuration q) {
  // distance du centre de c a q
  double distance = point_distance (c.xc, c.yc, q.x, q.y);
  //
  if ( fabs (distance - c.radius) > get_epsilon () ) { return false; }
  // angle de la droite orientee joignant le centre de c a q
  double angle = atan2 (q.y - c.yc, q.x - c.xc);
  // amgle de la direction de la mu-tangente
  if ( c.left && c.forward ) { angle = angle + Half_Pi - c.mu; }
  if ( c.left && !c.forward ) { angle = angle + Half_Pi + c.mu; }
  if ( !c.left && c.forward ) { angle = angle - Half_Pi + c.mu; }
  if ( !c.left && !c.forward ) { angle = angle - Half_Pi - c.mu; }
  angle = twopify (angle);
  //
  return fabs (q.theta - angle) < get_epsilon ();
}


// ############################################################################
// Class SMT_SC_Path
// ############################################################################
// constructeur
SMT_SC_Path::SMT_SC_Path (SMT_Configuration _start,
			  double _sharpness, bool _forward, double _length)
  : SMT_Path (_length, 0) {
  start = _start; sharpness = _sharpness; forward = _forward; }

// affichage alphanumerique des attributs
void SMT_SC_Path::print (bool eol) {
  cout << "* SC_Path: ";
  if ( forward ) cout << "forward"; else cout << "backward"; cout << " from ";
  start.print (false); 
  cout << ", sharpness " << sharpness;
  cout << ", length " << length;
  // eventuel saut de ligne final
  if ( eol ) { cout << endl; }
}


// ############################################################################
// Class SMT_CC_TP_Path
// ############################################################################
// constructeur
SMT_CC_TP_Path::SMT_CC_TP_Path (SMT_Configuration _start, 
				SMT_Configuration _end,
				double _kappa, double _sigma,
				SMT_Configuration *_qi1, 
				SMT_Configuration *_qi2,
				SMT_SC_Path *_sci1, 
				SMT_SC_Path *_sci2,
				bool forward, double _length, double _time) 
  : SMT_Path (_length, _time) {
  start = _start; end = _end; kappa = _kappa; sigma = _sigma;
  qi1 = _qi1; qi2 = _qi2; sci1 = _sci1; sci2 = _sci2; }

// destructeur
SMT_CC_TP_Path::~SMT_CC_TP_Path () {
  delete qi1; delete qi2; delete sci1; delete sci2; }

// affichage alphanumerique des attributs
void SMT_CC_TP_Path::print (bool eol) {
  cout << "* CC_TP_Path: ";
  start.print (false); cout << " -> "; 
  if ( !null (qi1) ) { qi1->print (false); cout << " -> "; }
  if ( !null (qi2) ) { qi2->print (false); cout << " -> "; }
  end.print (false); 
  cout << ", length " << length << ", time " << time; 
  if ( !null (sci1) ) sci1->print (false);
  if ( !null (sci2) ) sci2->print (false);
  // eventuel saut de ligne final
  if ( eol ) { cout << endl; }
}


