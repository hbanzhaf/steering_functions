// *** File: graphicize.cc
// *** Author(s): Th. Fraichard
// *** Last modified on 9 Jul 2008 

#include <graphicize.h>

// taille par defaut de la demi-branche de la croix representant
// graphiquement une configuration (repere "monde")
double Cross_Size = 0.05;

// gestion de la taille par defaut de la demi-branche de la croix
// representant graphiquement une configuration (repere "monde")
void set_cross_size (double size) { Cross_Size = size; }

// pas de discretisation par defaut pour l'affichage graphique d'une
// clotoide (repere "monde")
double Clothoid_Step = 0.05;

// gestion du pas de discretisation par defaut pour l'affichage graphique
// d'une clotoide (repere "monde")
void set_clothoid_step (double step) { Clothoid_Step = step; }

// transformation du cc-turn defini par un cc-circle c et une configuration
// d'arrivee q (supposee appartenir au cc-circle et avoir une orientation
// correcte) en un ensemble d'objets graphiques.  le booleen cusps indique
// si l'on se place dans le cas avec ou sans point de rebroussement
void cc_turn_graphicize (SMT_CC_Circle c, SMT_Configuration q, bool cusps,
			 FLTK_Display *display,
			 Fl_Color color, int style,  int width) {
  double deflection_min = twopify (pow (c.kappa, 2) / fabs (c.sigma));
  // calcul de la deflection
  double deflection;
  double alpha_c = c.start.theta;
  double alpha_q = q.theta;
  if ( c.left && c.forward ) { deflection = twopify (alpha_q - alpha_c); }
  if ( c.left && !c.forward ) { deflection = twopify (alpha_c - alpha_q); }
  if ( !c.left && c.forward ) { deflection = twopify (alpha_c - alpha_q); }
  if ( !c.left && !c.forward ) { deflection = twopify (alpha_q - alpha_c); }

  // *** deflection nulle -> cc-turn = segment
  FLTK_Line *li;
  if ( fabs (deflection) < get_epsilon () ) { 
    li = new FLTK_Line (c.start.x, c.start.y,
		       q.x, q.y,
		       color, style, width);
    display->add_object (li);
    return; }

  // *** cc-turn irregulier
  FLTK_Clothoid *cl;
  double x_i, y_i, theta_i, kappa_i;
  if ( deflection < deflection_min ) { 
    // sharpness et longueur d'un arc de clotoide
    double d1 = D1 (deflection / 2);
    double d2 = point_distance (c.start.x, c.start.y, 
				q.x, q.y);
    double sharpness = 4 * Pi * pow (d1, 2) / pow (d2, 2);
    double length =  sqrt (deflection / sharpness);
    // configuration intermediaire
    if ( !c.left ) { sharpness = - sharpness; } // signe de sigma
    // affichage arcs de clotoide
    end_of_clothoid (c.start.x, c.start.y,
		     c.start.theta, 0,
		     sharpness, c.forward, length,
		     &x_i, &y_i, &theta_i, &kappa_i);
    cl = new FLTK_Clothoid (c.start.x, c.start.y,
			   c.start.theta, 0,
			   x_i, y_i, theta_i, kappa_i,
			   sharpness, c.forward, length, Clothoid_Step,
			   color, style, width);
    display->add_object (cl);
    cl = new FLTK_Clothoid (q.x, q.y,
			   q.theta, 0,
			   x_i, y_i, theta_i, kappa_i,
			   sharpness, !c.forward, length, Clothoid_Step,
			   color, style, width);
    display->add_object (cl);
    return; }

  // *** cc-turn regulier
  FLTK_Arc *ar;
  double length_min = fabs (c.kappa / c.sigma);
  // affichage arcs de clotoide
  cl = new FLTK_Clothoid (c.start.x, c.start.y,
			 c.start.theta, 0,
			 c.first.x, c.first.y,
			 c.first.theta, c.first.kappa,
			 c.sigma, c.forward, length_min, Clothoid_Step,
			 color, style, width);
  display->add_object (cl);
  end_of_clothoid (q.x, q.y, q.theta, 0,
		   c.sigma, !c.forward, length_min,
		   &x_i, &y_i, &theta_i, &kappa_i);
  cl = new FLTK_Clothoid (q.x, q.y, q.theta, 0,
			 x_i, y_i, theta_i, kappa_i,
			 c.sigma, !c.forward, length_min, Clothoid_Step,
			 color, style, width);
  display->add_object (cl);
  // affichage rayons du cc-circle
  li = new FLTK_Line (c.xc, c.yc, c.first.x, c.first.y,
		     FL_GRAY, FL_DASH, 0);		
  display->add_object (li);
  li = new FLTK_Line (c.xc, c.yc, x_i, y_i, 
		     FL_GRAY, FL_DOT, 0);
  display->add_object (li);

  // cc-turn regulier
  bool ccw;
  if ( cusps && (deflection > deflection_min + Pi) ) {
    // cas avec manoeuvres
    if ( c.left ) { ccw = !c.forward; } else { ccw = c.forward; }
  }
  else {
    // cas sans manoeuvres
    if ( c.left ) { ccw = c.forward; } else { ccw = !c.forward; }
  }
  ar = new FLTK_Arc (c.first.x, c.first.y, x_i, y_i,
		    c.xc, c.yc, fabs (1 / c.kappa), ccw,
		    color, style, width);
  display->add_object (ar);
}

// configuration
void graphicize (SMT_Configuration q,
		 FLTK_Display *display,  
		 Fl_Color color, int style, int width) {
  double size = 0.1;
  FLTK_Cross *c;
  c = new FLTK_Cross (q.x, q.y, size, color, style, width);
  display->add_object (c);
  FLTK_Line *l;
  l = new FLTK_Line (q.x, q.y, q.x + 2 * size * cos (q.theta), 
		    q.y + 2 * size * sin (q.theta),
		    color, style, width);
  display->add_object (l);
}

// cc-circle
void graphicize (SMT_CC_Circle c,
		 FLTK_Display *display,  
		 Fl_Color color, int style, int width) {
  graphicize (c.start, display, color, style, width);
  graphicize (c.first, display, color, style, width);
  FLTK_Circle *circle;
  circle = new FLTK_Circle (c.xc, c.yc, c.radius, color, style, width);
  display->add_object (circle);
  FLTK_Line *line;
  line = new FLTK_Line (c.xc, c.yc, c.start.x, c.start.y, 
		       color, style, width);
  display->add_object (line); 
  FLTK_Clothoid *clothoid;
  clothoid = new FLTK_Clothoid (c.start.x, c.start.y, 
			       c.start.theta, c.start.kappa,
			       c.first.x, c.first.y, 
			       c.first.theta, c.first.kappa,
			       c.sigma, c.forward, 
			       fabs (c.kappa / c.sigma),
			       Clothoid_Step,
			       color, style, width);
  display->add_object (clothoid); 
}

// cc-dubins path
void graphicize (SMT_CC_Dubins_Path p,
		 FLTK_Display *display,  
		 Fl_Color color, int style, int width) {
  FLTK_Line *li;

  // configurations de depart et d'arrivee
  graphicize (p.start, display, color, style, width);
  graphicize (p.end, display, color, style, width);

  // configurations intermediaires
  if ( !null (p.qi1) ) 
    { graphicize (*(p.qi1), display, color, style, width); }
  if ( !null (p.qi2) ) 
    { graphicize (*(p.qi2), display, color, style, width); }

  switch ( p.type ) {
  case EMPTY: break;
  case S:
    li = new FLTK_Line (p.start.x, p.start.y,
		       p.end.x, p.end.y,
		       color, style, width);
    display->add_object (li);
    break;
  case R: case L:
    cc_turn_graphicize (*(p.cstart), p.end, false, 
			display, color, style, width);
    break;
  case SL: case SR:
    li = new FLTK_Line (p.start.x, p.start.y,
		       p.qi1->x, p.qi1->y,
		       color, style, width);
    display->add_object (li);
    cc_turn_graphicize (*(p.cend), *(p.qi1), false, 
			display, color, style, width);
    break;
  case RS: case LS:
    cc_turn_graphicize (*(p.cstart), *(p.qi1), false, 
			display, color, style, width);
    li = new FLTK_Line (p.qi1->x, p.qi1->y,
		       p.end.x, p.end.y,
		       color, style, width);
    display->add_object (li);
    break;
  case LSL: case LSR: case RSL: case RSR: 
    cc_turn_graphicize (*(p.cstart), *(p.qi1), false, 
			display, color, style, width);
    li = new FLTK_Line (p.qi1->x, p.qi1->y,
		       p.qi2->x, p.qi2->y,
		       color, style, width);
    display->add_object (li);
    cc_turn_graphicize (*(p.cend), *(p.qi2), false, 
			display, color, style, width);
    break;
  case LR1L: case LR2L: case RL1R: case RL2R:
    cc_turn_graphicize (*(p.cstart), *(p.qi1), false, 
			display, color, style, width);
    cc_turn_graphicize (*(p.ci1), *(p.qi2), false, 
			display, color, style, width);
    cc_turn_graphicize (*(p.cend), *(p.qi2), false, 
			display, color, style, width);
    break;
  default: break; }
}

// cc-reeds-shepp path
void graphicize (SMT_CC_RS_Path p,
		 FLTK_Display *display,  
		 Fl_Color color, int style, int width) {
  FLTK_Line *li;

  // configurations de depart et d'arrivee
  graphicize (p.start, display, color, style, width);
  graphicize (p.end, display, color, style, width);

  // configurations intermediaires
  if ( !null (p.qi1) ) 
    { graphicize (*(p.qi1), display, color, style, width); }
  if ( !null (p.qi2) ) 
    { graphicize (*(p.qi2), display, color, style, width); }
  if ( !null (p.qi3) ) 
    { graphicize (*(p.qi3), display, color, style, width); }
  if ( !null (p.qi4) ) 
    { graphicize (*(p.qi4), display, color, style, width); }

  switch ( p.type ) {
  case Empty: break;
  case Segment:
    li = new FLTK_Line (p.start.x, p.start.y,
		       p.end.x, p.end.y,
		       color, style, width);
    display->add_object (li);
    break;
  case T:
    cc_turn_graphicize (*(p.cstart), p.end, true, 
			display, color, style, width);
    break;
  case TT: case TcT:
    cc_turn_graphicize (*(p.cstart), *(p.qi1), true, 
			display, color, style, width);
    cc_turn_graphicize (*(p.cend), *(p.qi1), true, 
			display, color, style, width);
    break;
  case TST: case TcST: case TScT: case TcScT:
    cc_turn_graphicize (*(p.cstart), *(p.qi1), true, 
			display, color, style, width);
    li = new FLTK_Line (p.qi1->x, p.qi1->y,
		       p.qi2->x, p.qi2->y,
		       color, style, width);
    display->add_object (li);
    cc_turn_graphicize (*(p.cend), *(p.qi2), true, 
			display, color, style, width);
    break;
  case TTT: case TcTT: case TTcT: case TcTcT:
    cc_turn_graphicize (*(p.cstart), *(p.qi1), true, 
			display, color, style, width);
    cc_turn_graphicize (*(p.ci1), *(p.qi2), true, 
			display, color, style, width);
    cc_turn_graphicize (*(p.cend), *(p.qi2), true, 
			display, color, style, width);
    break;
  default: break; }
}

