// *** File: cc-rs.cc
// *** Author(s): Th. Fraichard
// *** Last modified on 3 May 2001

// Description: 

#include <cc-rs.h>

// ##### utilitaires ##########################################################

// ############################################################################
// pour chacune des familles des chemins de type cc_reeds_shepp, deux
// fonctions: la premiere (XXX_exists) verifie les conditions d'existence
// d'un chemin de la famille consideree, la deuxieme (XXX_path) le calcule.
// XXX_exists prend en entree deux cc-circles supposes decoules de la meme
// paire (kappa, sigma).  XXX_path suppose les conditions d'existence
// verifiees.


// ##### TT ###################################################################
bool TT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left == c2->left ) { return false; }
  if ( c1->forward == c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  if ( fabs (distance - 2 * c1->radius) > get_epsilon () ) { return false; }
  return true;
}

double TT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		SMT_Configuration **q) {
  // calcul de la position intermediaire
  double x = (c1->xc + c2->xc) / 2;
  double y = (c1->yc + c2->yc) / 2;
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double angle = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // calcul de l'orientation intermediaire
  double theta;
  if ( c1->left )
    if ( c1->forward ) { theta = angle + Half_Pi - c1->mu; }
    else { theta = angle + Half_Pi + c1->mu; }
  else 
    if ( c1->forward ) { theta = angle - Half_Pi + c1->mu; }
    else { theta = angle - Half_Pi - c1->mu; }
  // creation de la configuration intermediaire
  *q = new SMT_Configuration (x, y, theta, 0);
  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q) + rs_turn_length (c2, *q);
}


// ##### TcT ##################################################################
bool TcT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left == c2->left ) { return false; }
  if ( c1->forward != c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  if ( fabs (distance - 2 * c1->radius * cos (c1->mu)) > get_epsilon () ) 
    { return false; }
  return true;
}

double TcT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		 SMT_Configuration **q) {
  // distance des centres
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  // parametres de position de la configuration intermediaire dans le
  // repere local au centre de c1
  double delta_x = 0.5 * distance;
  double delta_y = sqrt (fabs (pow (c1->radius, 2) - pow (delta_x, 2)));
  //
  double x, y, theta;
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double angle = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // calcul de la configuration intermediaire
  if ( c1->left )
    if ( c1->forward ) { 
      theta = angle + Half_Pi; 
      global_frame_change (c1->xc, c1->yc, angle, delta_x, delta_y, &x, &y);
    }
    else { 
      theta = angle + Half_Pi; 
      global_frame_change (c1->xc, c1->yc, angle, delta_x, -delta_y, &x, &y); }
  else 
    if ( c1->forward ) { 
      theta = angle - Half_Pi;  
      global_frame_change (c1->xc, c1->yc, angle, delta_x, -delta_y, &x, &y); }
    else { 
      theta = angle - Half_Pi;  
      global_frame_change (c1->xc, c1->yc, angle, delta_x, delta_y, &x, &y); }

  // creation de la configuration intermediaire
  *q = new SMT_Configuration (x, y, theta, 0);
  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q) + rs_turn_length (c2, *q);
}


// ##### TST ##################################################################
bool TiST_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left == c2->left ) { return false; }
  if ( c1->forward == c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return ( distance >= 2 * c1->radius); 
}   
    
bool TeST_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left != c2->left ) { return false; }
  if ( c1->forward == c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return ( distance >= 2 * c1->radius * sin (c1->mu)); 
}   
    
bool TST_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  return TiST_exists (c1, c2) || TeST_exists (c1, c2);
}    
     
double TiST_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		  SMT_Configuration **q1, SMT_Configuration **q2) {
  // distance des centres
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double angle = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // angle entre la droite des centres et la direction de la mu-tangente
  double alpha = fabs (asin (2 * c1->radius * cos (c1->mu) / distance));
  // parametres de position des configurations intermediaires dans le repere
  // local au centre de c1
  double delta_x = fabs (c1->radius * sin (c1->mu));
  double delta_y = fabs (c1->radius * cos (c1->mu));
  //
  double x, y, theta;
  if ( c1->left && c1->forward ) {
    theta = angle + alpha;
    global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( c1->left && !c1->forward ) {
    theta = angle - alpha;
    global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( !c1->left && c1->forward ) {
    theta = angle - alpha;
    global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( !c1->left && !c1->forward ) {
    theta = angle + alpha;
    global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }

  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q1) 
    + configuration_distance (**q1, **q2) 
    + rs_turn_length (c2, *q2);
}   
    
double TeST_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		  SMT_Configuration **q1, SMT_Configuration **q2) {
  // parametres de position des configurations intermediaires dans le
  // repere local au centre de c1
  double delta_x = fabs (c1->radius * sin (c1->mu));
  double delta_y = fabs (c1->radius * cos (c1->mu));
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  //
  double x, y;
  // calcul et creation des configurations intermediaires
  if ( c1->left && c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( c1->left && !c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( !c1->left && c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( !c1->left && !c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q1) 
    + configuration_distance (**q1, **q2) 
    + rs_turn_length (c2, *q2);
}   

double TST_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		 SMT_Configuration **q1, SMT_Configuration **q2) {
  if ( TiST_exists (c1, c2) ) { return TiST_path (c1, c2, q1, q2); }
  if ( TeST_exists (c1, c2) ) { return TeST_path (c1, c2, q1, q2); }
  return FLT_MAX;		// pour eviter un warning a la compilation
}   


// ##### TcST ################################################################
bool TciST_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left == c2->left ) { return false; }
  if ( c1->forward != c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance >= 2 * c1->radius * cos (c1->mu); 
}   
    
bool TceST_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left != c2->left ) { return false; }
  if ( c1->forward != c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance >= get_epsilon (); 
}   
    
bool TcST_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  return TciST_exists (c1, c2) || TceST_exists (c1, c2);
}    
     
double TciST_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		   SMT_Configuration **q1, SMT_Configuration **q2) {
  // distance des centres
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double angle = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // angle entre la droite des centres et la direction de la mu-tangente
  double alpha = fabs (asin (2 * c1->radius * cos (c1->mu) / distance));
  // parametres de position des configurations intermediaires dans le repere
  // local au centre de c1
  double delta_x = fabs (c1->radius * sin (c1->mu));
  double delta_y = fabs (c1->radius * cos (c1->mu));
  //
  double x, y, theta;
  if ( c1->left && c1->forward ) {
    theta = angle - alpha;
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( c1->left && !c1->forward ) {
    theta = angle + alpha;
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( !c1->left && c1->forward ) {
    theta = angle + alpha;
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( !c1->left && !c1->forward ) {
    theta = angle - alpha;
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }

  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q1) 
    + configuration_distance (**q1, **q2) 
    + rs_turn_length (c2, *q2);
}   
    
double TceST_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		  SMT_Configuration **q1, SMT_Configuration **q2) {
  // parametres de position des configurations intermediaires dans le
  // repere local au centre de c1
  double delta_x = fabs (c1->radius * sin (c1->mu));
  double delta_y = fabs (c1->radius * cos (c1->mu));
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  //
  double x, y;
  // calcul et creation des configurations intermediaires
  if ( c1->left && c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( c1->left && !c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( !c1->left && c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( !c1->left && !c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, -delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q1) 
    + configuration_distance (**q1, **q2) 
    + rs_turn_length (c2, *q2);
}   

double TcST_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		 SMT_Configuration **q1, SMT_Configuration **q2) {
  if ( TciST_exists (c1, c2) ) { return TciST_path (c1, c2, q1, q2); }
  if ( TceST_exists (c1, c2) ) { return TceST_path (c1, c2, q1, q2); }
  return FLT_MAX;		// pour eviter un warning a la compilation
}   


// ##### TScT #################################################################
bool TiScT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left == c2->left ) { return false; }
  if ( c1->forward != c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance >= 2 * c1->radius * cos (c1->mu); 
}
    
bool TeScT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left != c2->left ) { return false; }
  if ( c1->forward != c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance >= get_epsilon (); 
}   
    
bool TScT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  return TiScT_exists (c1, c2) || TeScT_exists (c1, c2);
}    
     
double TiScT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		   SMT_Configuration **q1, SMT_Configuration **q2) {
  // distance des centres
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double angle = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // angle entre la droite des centres et la direction de la mu-tangente
  double alpha = fabs (asin (2 * c1->radius * cos (c1->mu) / distance));
  // parametres de position des configurations intermediaires dans le repere
  // local au centre de c1
  double delta_x = fabs (c1->radius * sin (c1->mu));
  double delta_y = fabs (c1->radius * cos (c1->mu));
  //
  double x, y, theta;
  if ( c1->left && c1->forward ) {
    theta = angle + alpha;
    global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( c1->left && !c1->forward ) {
    theta = angle - alpha;
    global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( !c1->left && c1->forward ) {
    theta = angle - alpha;
    global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( !c1->left && !c1->forward ) {
    theta = angle + alpha;
    global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }

  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q1) 
    + configuration_distance (**q1, **q2) 
    + rs_turn_length (c2, *q2);
}   
    
double TeScT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		  SMT_Configuration **q1, SMT_Configuration **q2) {
  // parametres de position des configurations intermediaires dans le
  // repere local au centre de c1
  double delta_x = fabs (c1->radius * sin (c1->mu));
  double delta_y = fabs (c1->radius * cos (c1->mu));
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  //
  double x, y;
  // calcul et creation des configurations intermediaires
  if ( c1->left && c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( c1->left && !c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( !c1->left && c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( !c1->left && !c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q1) 
    + configuration_distance (**q1, **q2) 
    + rs_turn_length (c2, *q2);
}   

double TScT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		 SMT_Configuration **q1, SMT_Configuration **q2) {
  if ( TiScT_exists (c1, c2) ) { return TiScT_path (c1, c2, q1, q2); }
  if ( TeScT_exists (c1, c2) ) { return TeScT_path (c1, c2, q1, q2); }
  return FLT_MAX;		// pour eviter un warning a la compilation
}   


// ##### TcScT ################################################################
bool TciScT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left == c2->left ) { return false; }
  if ( c1->forward == c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance >= 2 * c1->radius * cos (c1->mu); 
}
    
bool TceScT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left != c2->left ) { return false; }
  if ( c1->forward == c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance >= get_epsilon (); 
}   
    
bool TcScT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  return TciScT_exists (c1, c2) || TceScT_exists (c1, c2);
}    
     
double TciScT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		   SMT_Configuration **q1, SMT_Configuration **q2) {
  // distance des centres
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double angle = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // angle entre la droite des centres et la direction de la mu-tangente
  double alpha = fabs (asin (2 * c1->radius * cos (c1->mu) / distance));
  // parametres de position des configurations intermediaires dans le repere
  // local au centre de c1
  double delta_x = fabs (c1->radius * sin (c1->mu));
  double delta_y = fabs (c1->radius * cos (c1->mu));
  //
  double x, y, theta;
  if ( c1->left && c1->forward ) {
    theta = angle - alpha;
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( c1->left && !c1->forward ) {
    theta = angle + alpha;
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( !c1->left && c1->forward ) {
    theta = angle + alpha;
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( !c1->left && !c1->forward ) {
    theta = angle - alpha;
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }

  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q1) 
    + configuration_distance (**q1, **q2) 
    + rs_turn_length (c2, *q2);
}   
    
double TceScT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		    SMT_Configuration **q1, SMT_Configuration **q2) {
  // parametres de position des configurations intermediaires dans le
  // repere local au centre de c1
  double delta_x = fabs (c1->radius * sin (c1->mu));
  double delta_y = fabs (c1->radius * cos (c1->mu));
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  //
  double x, y;
  // calcul et creation des configurations intermediaires
  if ( c1->left && c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( c1->left && !c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  if ( !c1->left && c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, -delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta + Pi, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, -delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta + Pi, 0);
  }
  if ( !c1->left && !c1->forward ) {
    global_frame_change (c1->xc, c1->yc, theta, -delta_x, delta_y, &x, &y);
    *q1 = new SMT_Configuration (x, y, theta, 0);
    global_frame_change (c2->xc, c2->yc, theta, delta_x, delta_y, &x, &y);
    *q2 = new SMT_Configuration (x, y, theta, 0);
  }
  // calcul de la longueur du chemin
  return rs_turn_length (c1, *q1) 
    + configuration_distance (**q1, **q2) 
    + rs_turn_length (c2, *q2);
}   

double TcScT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		   SMT_Configuration **q1, SMT_Configuration **q2) {
  if ( TciScT_exists (c1, c2) ) { return TciScT_path (c1, c2, q1, q2); }
  if ( TceScT_exists (c1, c2) ) { return TceScT_path (c1, c2, q1, q2); }
  return FLT_MAX;		// pour eviter un warning a la compilation
}   


// ##### TTT ##################################################################
void TTT_tangent_circles (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
			  SMT_Configuration **q1, SMT_Configuration **q2,
			  SMT_Configuration **q3, SMT_Configuration **q4) {
  // distance des centres
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // parametres de position des centres des cc-circles tangents dans le
  // repere local au centre de c1
  double r = 2 * c1->radius;
  double delta_x = 0.5 * distance;
  double delta_y = sqrt (fabs (pow (delta_x, 2) - pow (r, 2)));

  // creation des cc-circles tangents
  double x, y;
  global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
  SMT_CC_Circle tgt1 (x, y, c1->radius, !c1->left, c1->forward,
		      c1->kappa, c1->sigma, c1->mu);
  global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
  SMT_CC_Circle tgt2 (x, y, c1->radius, !c1->left, c1->forward,
		      c1->kappa, c1->sigma, c1->mu);
  // calcul des configurations intermediaires
  TT_path (c1, &tgt1, q1);
  TT_path (&tgt1, c2, q2);
  TT_path (c1, &tgt2, q3);
  TT_path (&tgt2, c2, q4);
}

bool TTT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left != c2->left ) { return false; }
  if ( c1->forward == c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance <= 4 * c1->radius;   
}
     
double TTT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		 SMT_Configuration **q1, SMT_Configuration **q2,
		 SMT_CC_Circle **ci) {
  // 
  SMT_Configuration *qa, *qb, *qc, *qd;
  TTT_tangent_circles (c1, c2, &qa, &qb, &qc, &qd);
  //
  SMT_CC_Circle *middle1, *middle2;
  middle1 = new SMT_CC_Circle (*qa, !c1->left, c1->forward, 
			       c1->kappa, c1->sigma);
  double length1 = rs_turn_length (c1, qa) 
    + rs_turn_length (middle1, qb) 
    + rs_turn_length (c2, qb);
  middle2 = new SMT_CC_Circle (*qc, !c1->left, c1->forward, 
			       c1->kappa, c1->sigma);
  double length2 = rs_turn_length (c1, qc) 
    + rs_turn_length (middle2, qd) 
    + rs_turn_length (c2, qd);
  if ( length1 < length2 ) {
    *q1 = qa; *q2 = qb; *ci = middle1; 
    delete qc; delete qd; delete middle2; 
    return length1;
  }
  else {
    *q1 = qc; *q2 = qd; *ci = middle2; 
    delete qa; delete qb; delete middle1; 
    return length2;
  }
  return FLT_MAX;		// pour eviter un warning a la compilation
}


// ##### TcTT #################################################################
void TcTT_tangent_circles (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
			   SMT_Configuration **q1, SMT_Configuration **q2,
			   SMT_Configuration **q3, SMT_Configuration **q4) {
  // distance des centres
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // parametres de position des centres des cc-circles tangents dans le
  // repere local au centre de c1
  double r1 = 2 * c1->radius * cos (c1->mu);
  double r2 = 2 * c1->radius;
  double delta_x = (pow (r1, 2) + pow (distance, 2) - pow (r2, 2)) 
    / (2 * distance);
  double delta_y = sqrt (fabs (pow (r1, 2) - pow (delta_x, 2)));

  // creation des cc-circles tangents
  double x, y;
  global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
  SMT_CC_Circle tgt1 (x, y, c1->radius, !c1->left, !c1->forward,
		      c1->kappa, c1->sigma, c1->mu);
  global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
  SMT_CC_Circle tgt2 (x, y, c1->radius, !c1->left, !c1->forward,
		      c1->kappa, c1->sigma, c1->mu);
  // calcul des configurations intermediaires
  TcT_path (c1, &tgt1, q1);
  TT_path (&tgt1, c2, q2);
  TcT_path (c1, &tgt2, q3);
  TT_path (&tgt2, c2, q4);
}

bool TcTT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left != c2->left ) { return false; }
  if ( c1->forward != c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance <= 2 * c1->radius * (1 + cos (c1->mu));   
}
     
double TcTT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		  SMT_Configuration **q1, SMT_Configuration **q2,
		  SMT_CC_Circle **ci) {
  // 
  SMT_Configuration *qa, *qb, *qc, *qd;
  TcTT_tangent_circles (c1, c2, &qa, &qb, &qc, &qd);
  //
  SMT_CC_Circle *middle1, *middle2;
  middle1 = new SMT_CC_Circle (*qa, !c1->left, !c1->forward, 
			       c1->kappa, c1->sigma);
  double length1 = rs_turn_length (c1, qa) 
    + rs_turn_length (middle1, qb) 
    + rs_turn_length (c2, qb);
  middle2 = new SMT_CC_Circle (*qc, !c1->left, !c1->forward, 
			       c1->kappa, c1->sigma);
  double length2 = rs_turn_length (c1, qc) 
    + rs_turn_length (middle2, qd) 
    + rs_turn_length (c2, qd);
  if ( length1 < length2 ) {
    *q1 = qa; *q2 = qb; *ci = middle1; 
    delete qc; delete qd; delete middle2;
    return length1;
  }
  else {
    *q1 = qc; *q2 = qd; *ci = middle2; 
    delete qa; delete qb; delete middle1; 
    return length2;
  }
  return FLT_MAX;		// pour eviter un warning a la compilation
}


// ##### TTcT #################################################################
void TTcT_tangent_circles (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
			   SMT_Configuration **q1, SMT_Configuration **q2,
			   SMT_Configuration **q3, SMT_Configuration **q4) {
  // distance des centres
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // parametres de position des centres des cc-circles tangents dans le
  // repere local au centre de c1
  double r1 = 2 * c1->radius;
  double r2 = 2 * c1->radius * cos (c1->mu);
  double delta_x = (pow (r1, 2) + pow (distance, 2) - pow (r2, 2)) 
    / (2 * distance);
  double delta_y = sqrt (fabs (pow (r1, 2) - pow (delta_x, 2)));

  // creation des cc-circles tangents
  double x, y;
  global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
  SMT_CC_Circle tgt1 (x, y, c1->radius, !c1->left, c1->forward,
		      c1->kappa, c1->sigma, c1->mu);
  global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
  SMT_CC_Circle tgt2 (x, y, c1->radius, !c1->left, c1->forward,
		      c1->kappa, c1->sigma, c1->mu);
  // calcul des configurations intermediaires
  TT_path (c1, &tgt1, q1);
  TcT_path (&tgt1, c2, q2);
  TT_path (c1, &tgt2, q3);
  TcT_path (&tgt2, c2, q4);
}

bool TTcT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left != c2->left ) { return false; }
  if ( c1->forward != c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance <= 2 * c1->radius * (1 + cos (c1->mu));   
}
     
double TTcT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		 SMT_Configuration **q1, SMT_Configuration **q2,
		 SMT_CC_Circle **ci) {
  // 
  SMT_Configuration *qa, *qb, *qc, *qd;
  TTcT_tangent_circles (c1, c2, &qa, &qb, &qc, &qd);
  //
  SMT_CC_Circle *middle1, *middle2;
  middle1 = new SMT_CC_Circle (*qa, !c1->left, c1->forward, 
			       c1->kappa, c1->sigma);
  double length1 = rs_turn_length (c1, qa) 
    + rs_turn_length (middle1, qb) 
    + rs_turn_length (c2, qb);
  middle2 = new SMT_CC_Circle (*qc, !c1->left, c1->forward, 
			       c1->kappa, c1->sigma);
  double length2 = rs_turn_length (c1, qc) 
    + rs_turn_length (middle2, qd) 
    + rs_turn_length (c2, qd);
  if ( length1 < length2 ) {
    *q1 = qa; *q2 = qb; *ci = middle1; 
    delete qc; delete qd; delete middle2; 
    return length1;
  }
  else {
    *q1 = qc; *q2 = qd; *ci = middle2; 
    delete qa; delete qb; delete middle1; 
    return length2;
  }
  return FLT_MAX;		// pour eviter un warning a la compilation
}


// ##### TcTcT ################################################################
void TcTcT_tangent_circles (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
			    SMT_Configuration **q1, SMT_Configuration **q2,
			    SMT_Configuration **q3, SMT_Configuration **q4) {
  // distance des centres
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  // angle de la droite orientee joignant le centre de c1 au centre de c2
  double theta = atan2 (c2->yc - c1->yc, c2->xc - c1->xc);
  // parametres de position des centres des cc-circles tangents dans le
  // repere local au centre de c1
  double r = 2 * c1->radius * cos (c1->mu);
  double delta_x = 0.5 * distance;
  double delta_y = sqrt (fabs (pow (r, 2) - pow (delta_x, 2)));

  // creation des cc-circles tangents
  double x, y;
  global_frame_change (c1->xc, c1->yc, theta, delta_x, delta_y, &x, &y);
  SMT_CC_Circle tgt1 (x, y, c1->radius, !c1->left, !c1->forward,
		      c1->kappa, c1->sigma, c1->mu);
  global_frame_change (c1->xc, c1->yc, theta, delta_x, -delta_y, &x, &y);
  SMT_CC_Circle tgt2 (x, y, c1->radius, !c1->left, !c1->forward,
		      c1->kappa, c1->sigma, c1->mu);
  // calcul des configurations intermediaires
  TcT_path (c1, &tgt1, q1);
  TcT_path (&tgt1, c2, q2);
  TcT_path (c1, &tgt2, q3);
  TcT_path (&tgt2, c2, q4);
}

bool TcTcT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left != c2->left ) { return false; }
  if ( c1->forward == c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  return distance <= 4 * c1->radius * cos (c1->mu);   
}
     
double TcTcT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		   SMT_Configuration **q1, SMT_Configuration **q2,
		   SMT_CC_Circle **ci) {
  // 
  SMT_Configuration *qa, *qb, *qc, *qd;
  TcTcT_tangent_circles (c1, c2, &qa, &qb, &qc, &qd);
  //
  SMT_CC_Circle *middle1, *middle2;
  middle1 = new SMT_CC_Circle (*qa, !c1->left, !c1->forward, 
			       c1->kappa, c1->sigma);
  double length1 = rs_turn_length (c1, qa) 
    + rs_turn_length (middle1, qb) 
    + rs_turn_length (c2, qb);
  middle2 = new SMT_CC_Circle (*qc, !c1->left, !c1->forward, 
			       c1->kappa, c1->sigma);
  double length2 = rs_turn_length (c1, qc) 
    + rs_turn_length (middle2, qd) 
    + rs_turn_length (c2, qd);
  if ( length1 < length2 ) {
    *q1 = qa; *q2 = qb; *ci = middle1; 
    delete qc; delete qd; delete middle2; 
    return length1;
  }
  else {
    *q1 = qc; *q2 = qd; *ci = middle2; 
    delete qa; delete qb; delete middle1; 
    return length2;
  }
  return FLT_MAX;		// pour eviter un warning a la compilation
}


// ##### TTST #################################################################
// ##### TcTST ################################################################
// ##### TTcST ################################################################
// ##### TTScT ################################################################
// ##### TcTcST ###############################################################
// ##### TcTScT ###############################################################
// ##### TTcScT ###############################################################
// ##### TcTcScT ##############################################################


// ##### TTSTT ################################################################
// ##### TcTSTT ###############################################################
// ##### TTcSTT ###############################################################
// ##### TTScTT ###############################################################
// ##### TTSTcT ###############################################################
// ##### TcTcSTT ##############################################################
// ##### TTcScTT ##############################################################
// ##### TTScTcT ##############################################################
// ##### TcTScTT ##############################################################
// ##### TcTSTcT ##############################################################
// ##### TTcSTcT ##############################################################
// ##### TcTcScTT #############################################################
// ##### TcTcSTcT #############################################################
// ##### TcTScTcT #############################################################
// ##### TTcScTcT #############################################################
// ##### TcTcScTcT ############################################################


// ##### TTSTT ################################################################
bool TTiSTT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left == c2->left ) { return false; }
  if ( c1->forward == c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  if ( fabs (distance - 6 * c1->radius) > get_epsilon () ) { return false; }
  return true;
}

bool TTeSTT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  // les cc-circles sont-ils compatibles ?
  if ( c1->left != c2->left ) { return false; }
  if ( c1->forward == c2->forward ) { return false; }
  // les centres sont-ils a la bonne distance ?
  double distance = point_distance (c1->xc, c1->yc, c2->xc, c2->yc);
  if ( fabs (distance - 2 * c1->radius * (2 + sin (c1->mu))) 
       > get_epsilon () ) 
    { return false; }
  return true;
}

bool TTSTT_exists (SMT_CC_Circle *c1, SMT_CC_Circle *c2) {
  return TTiSTT_exists (c1, c2) || TTeSTT_exists (c1, c2);
}

double TTSTT_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
		   SMT_Configuration **q1, SMT_Configuration **q2,
		   SMT_Configuration **q3, SMT_Configuration **q4,
		   SMT_CC_Circle **ci1, SMT_CC_Circle **ci2) {
  return FLT_MAX;		// pour eviter un warning a la compilation
}

// ############################################################################
SMT_CC_RS_Path *cc_circles_rs_path (SMT_CC_Circle *c1, SMT_CC_Circle *c2,
				    double kappa, double sigma) {
  // tableaux destines a stocker les longueurs des chemins possibles, les
  // configurations intermediaires et les cc-circles correspondants
  double length [nb_cc_rs_paths];
  double_array_init (length, nb_cc_rs_paths, FLT_MAX);
  SMT_Configuration *qi1 [nb_cc_rs_paths];
  pointer_array_init ((void **) qi1, nb_cc_rs_paths);
  SMT_Configuration *qi2 [nb_cc_rs_paths];
  pointer_array_init ((void **) qi2, nb_cc_rs_paths);
  SMT_Configuration *qi3 [nb_cc_rs_paths];
  pointer_array_init ((void **) qi3, nb_cc_rs_paths);
  SMT_Configuration *qi4 [nb_cc_rs_paths];
  pointer_array_init ((void **) qi4, nb_cc_rs_paths);
  SMT_CC_Circle *cstart [nb_cc_rs_paths];
  pointer_array_init ((void **) cstart, nb_cc_rs_paths);
  SMT_CC_Circle *ci1 [nb_cc_rs_paths];
  pointer_array_init ((void **) ci1, nb_cc_rs_paths);
  SMT_CC_Circle *ci2 [nb_cc_rs_paths];
  pointer_array_init ((void **) ci2, nb_cc_rs_paths);
  SMT_CC_Circle *cend [nb_cc_rs_paths];
  pointer_array_init ((void **) cend, nb_cc_rs_paths);

  // traitement des differents cas possibles
  // d'abord les cas particuliers
  // cas Empty
  if ( configuration_equal (c1->start, c2->start) ) {
    length [Empty] = 0;
    goto label_fin;
  }
  // cas Segment+
  if (configuration_aligned (c1->start, c2->start) ) {
    length [Segment] = configuration_distance (c1->start, c2->start);
    goto label_fin;
  }
  // cas Segment-
  if (configuration_aligned (c2->start, c1->start) ) {
    length [Segment] = configuration_distance (c2->start, c1->start);
    goto label_fin;
  }
  // cas T
  if ( configuration_on_cc_circle (*c1, c2->start) ) {
    cstart [T] = new SMT_CC_Circle (*c1);
    length [T] = rs_turn_length (c1, &(c2->start));
    goto label_fin;
  }

  // puis les autres
  // cas TT
  // ######
  if ( TT_exists (c1, c2) ) { 
    cstart [TT] = new SMT_CC_Circle (*c1);
    cend [TT] = new SMT_CC_Circle (*c2);
    length [TT] = TT_path (c1, c2, &qi1 [TT]);
  }
  // cas TcT
  // #######
  if ( TcT_exists (c1, c2) ) {
    cstart [TcT] = new SMT_CC_Circle (*c1);
    cend [TcT] = new SMT_CC_Circle (*c2);
    length [TcT] = TcT_path (c1, c2, &qi1 [TcT]);
  }

  // cas TST
  // #######
  if ( TST_exists (c1, c2) ) {
    cstart [TST] = new SMT_CC_Circle (*c1);
    cend [TST] = new SMT_CC_Circle (*c2);
    length [TST] = TST_path (c1, c2, &qi1 [TST], &qi2 [TST]);
  }
  // cas TcST
  // ########
  if ( TcST_exists (c1, c2) ) {
    cstart [TcST] = new SMT_CC_Circle (*c1);
    cend [TcST] = new SMT_CC_Circle (*c2);
    length [TcST] = TcST_path (c1, c2, &qi1 [TcST], &qi2 [TcST]);
  }
  // cas TScT
  // ########
  if ( TScT_exists (c1, c2) ) {
    cstart [TScT] = new SMT_CC_Circle (*c1);
    cend [TScT] = new SMT_CC_Circle (*c2);
    length [TScT] = TScT_path (c1, c2, &qi1 [TScT], &qi2 [TScT]);
  }
  // cas TcScT
  // #########
  if ( TcScT_exists (c1, c2) ) {
    cstart [TcScT] = new SMT_CC_Circle (*c1);
    cend [TcScT] = new SMT_CC_Circle (*c2);
    length [TcScT] = TcScT_path (c1, c2, &qi1 [TcScT], &qi2 [TcScT]);
  }

  // cas TTT
  // #######
  if ( TTT_exists (c1, c2) ) {
    cstart [TTT] = new SMT_CC_Circle (*c1);
    cend [TTT] = new SMT_CC_Circle (*c2);
    length [TTT] = TTT_path (c1, c2, &qi1 [TTT], &qi2 [TTT], &ci1 [TTT]);
  }
  // cas TcTT
  // ########
  if ( TcTT_exists (c1, c2) ) {
    cstart [TcTT] = new SMT_CC_Circle (*c1);
    cend [TcTT] = new SMT_CC_Circle (*c2);
    length [TcTT] = TcTT_path (c1, c2, &qi1 [TcTT], &qi2 [TcTT], &ci1 [TcTT]);
  }
  // cas TTcT
  // ########
  if ( TTcT_exists (c1, c2) ) {
    cstart [TTcT] = new SMT_CC_Circle (*c1);
    cend [TTcT] = new SMT_CC_Circle (*c2);
    length [TTcT] = TTcT_path (c1, c2, &qi1 [TTcT], &qi2 [TTcT], &ci1 [TTcT]);
  }
  // cas TcTcT
  // #########
  if ( TcTcT_exists (c1, c2) ) {
    cstart [TcTcT] = new SMT_CC_Circle (*c1);
    cend [TcTcT] = new SMT_CC_Circle (*c2);
    length [TcTcT] = TcTcT_path (c1, c2, &qi1 [TcTcT], 
				 &qi2 [TcTcT], &ci1 [TcTcT]);
  }

 label_fin:  
  // determination du plus court chemin
  cc_rs_path_type best_path = (cc_rs_path_type) 
    array_index_min (length, nb_cc_rs_paths);

  // contruction du plus court chemin
  SMT_CC_RS_Path *path;
  path = new SMT_CC_RS_Path (c1->start, c2->start,
			     best_path, kappa, sigma,
			     qi1 [best_path], qi2 [best_path],
			     qi3 [best_path], qi4 [best_path],
			     cstart [best_path], cend [best_path],
			     ci1 [best_path], ci2 [best_path],
			     length [best_path], 0);

  // menage
  for (int i = 0; i < nb_cc_rs_paths; i++) {
    if ( i != best_path ) {
      delete qi1 [i]; delete qi2 [i]; delete qi3 [i]; delete qi4 [i];
      delete cstart [i]; delete ci1 [i]; delete ci2 [i]; delete cend [i];
    }
  }

  //retour du plus court chemin
  return path;
}

// ############################################################################
SMT_CC_RS_Path *cc_reeds_shepp (SMT_Configuration start, 
				SMT_Configuration end,
				double kappa, double sigma) {
  // on lance le chronometre
  start_chronometer ();

  // calcul des cc_circles associes aux configurations start et end
  SMT_CC_Circle *start_circle [4];
  SMT_CC_Circle *end_circle [4];
  
  start_circle [0] = new SMT_CC_Circle (start, true, true, kappa, sigma);
  start_circle [1] = new SMT_CC_Circle (start, false, true, kappa, sigma);
  start_circle [2] = new SMT_CC_Circle (start, true, false, kappa, sigma);
  start_circle [3] = new SMT_CC_Circle (start, false, false, kappa, sigma);
  end_circle [0] = new SMT_CC_Circle (end, true, true, kappa, sigma);
  end_circle [1] = new SMT_CC_Circle (end, false, true, kappa, sigma);
  end_circle [2] = new SMT_CC_Circle (end, true, false, kappa, sigma);
  end_circle [3] = new SMT_CC_Circle (end, false, false, kappa, sigma);

  // calcul du cc-rs-path pour les 16 combinaisons possibles de cc-circles
  // de debut et de fin
  SMT_CC_RS_Path *path [] = {0, 0, 0, 0, 0, 0, 0, 0, 
			     0, 0, 0, 0, 0, 0, 0, 0};
  double lg [] = {FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX,
		  FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX,
		  FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX,
		  FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX};
    
  for (int i = 0; i < 4; i++) {	// index des cc-circles de depart
    for (int j = 0; j < 4; j++) { // index des cc-circles d'arrivee
      path [4*i + j] = cc_circles_rs_path (start_circle [i], end_circle [j],
					   kappa, sigma);
      if ( !null (path [4*i + j]) ) {lg [4*i + j] = 
				       path [4*i + j]->length; }
    }
  }

  // determination du plus court chemin
  int best_path = array_index_min (lg, 15);
  
  // mesure du temps ecoule
  double time = look_at_chronometer ();
  path [best_path]->time = time;

  // trace des calculs
//    cout << "*** CC_Reeds_Shepp" << endl;
//    for (int i = 0; i < 16; i++) {
//      cout << i << ": "; 
//      if ( !null (path [i]) ) { path [i]->print (true); }
//    }  
//    cout << "plus court chemin: " << (int) best_path << endl;

  // menage
  for (int i = 0; i < 4; i++) { 
    delete start_circle [i]; delete end_circle [i];
  }
  for (int i = 0; i < 16; i++) { 
    if ( i != best_path) { delete path [i]; }
  }

  //retour du plus court chemin
  return path [best_path];
}
