// *** File: cc-dubins.cc
// *** Author(s): Th. Fraichard
// *** Last modified on 5 Sep 2000

// Description: 

#include <cc-dubins.h>

// note: faire deux fois delete sur le meme objet a des consequences
// desastreuses.  ceci explique les creations multiples d'un meme objet.
// c'est un peu lourd mais cela permet d'eviter les problemes au moment du
// menage...

SMT_CC_Dubins_Path *cc_dubins (SMT_Configuration start, 
			       SMT_Configuration end,
			       double kappa, double sigma) {
  // on lance le chronometre
  start_chronometer ();

  // tableaux destines a stocker les longueurs des chemins possibles, les
  // configurations intermediaires et les cc-circles correspondants
  double length [nb_cc_dubins_paths];
  double_array_init (length, nb_cc_dubins_paths, FLT_MAX);
  SMT_Configuration *qi1 [nb_cc_dubins_paths];
  pointer_array_init ((void **) qi1, nb_cc_dubins_paths);
  SMT_Configuration *qi2 [nb_cc_dubins_paths];
  pointer_array_init ((void **) qi2, nb_cc_dubins_paths);
  SMT_CC_Circle *cstart [nb_cc_dubins_paths];
  pointer_array_init ((void **) cstart, nb_cc_dubins_paths);
  SMT_CC_Circle *ci1 [nb_cc_dubins_paths];
  pointer_array_init ((void **) ci1, nb_cc_dubins_paths);
  SMT_CC_Circle *cend [nb_cc_dubins_paths];
  pointer_array_init ((void **) cend, nb_cc_dubins_paths);

  // calcul des cc_circles associes aux configurations start et end
  SMT_CC_Circle *start_left_forward, *start_right_forward;
  SMT_CC_Circle *end_left_backward, *end_right_backward; 
  
  start_left_forward = new SMT_CC_Circle (start, true, true, kappa, sigma);
  start_right_forward = new SMT_CC_Circle (start, false, true, kappa, sigma);
  end_left_backward = new SMT_CC_Circle (end, true, false, kappa, sigma);
  end_right_backward = new SMT_CC_Circle (end, false, false, kappa, sigma);

  // traitement des differents cas possibles
  // d'abord les cas particuliers
  // cas EMPTY
  if ( configuration_equal (start, end) ) {
    length [EMPTY] = 0;
    goto label_fin;
  }
  // cas S
  if (configuration_aligned (start, end) ) {
    length [S] = configuration_distance (start, end);
    goto label_fin;
  }
  // cas L
  if ( configuration_on_cc_circle (*start_left_forward, end) ) {
    cstart [L] = new SMT_CC_Circle (*start_left_forward);
    length [L] = dubins_turn_length (start_left_forward, &end);
    goto label_fin;
  }
  // cas R
  if ( configuration_on_cc_circle (*start_right_forward, end) ) {
    cstart [R] = new SMT_CC_Circle (*start_right_forward);
    length [R] = dubins_turn_length (start_right_forward, &end);
    goto label_fin;
  }

  // puis les autres
  // cas LSL et sous-cas LS, SL
  SMT_Configuration *qa, *qb;
  if ( external_mu_tangent_exists (*start_left_forward, *end_left_backward) )
    {
      external_mu_tangent (*start_left_forward, *end_left_backward, 
			   &qa, &qb);
      // sous-cas LS
      if ( configuration_aligned (*qb, end) ) {
	cstart [LS] = new SMT_CC_Circle (*start_left_forward);
	qi1 [LS] = new SMT_Configuration (*qa);
	length [LS] = dubins_turn_length (start_left_forward, qa)
	  + configuration_distance (*qa, end); }
      // sous-cas SL
      else if ( configuration_aligned (start, *qa) ) {
	cend [SL] = new SMT_CC_Circle (*end_left_backward);
	qi1 [SL] = new SMT_Configuration (*qb);
	length [SL] = configuration_distance (start, *qb)
	  + dubins_turn_length (end_left_backward, qb); }
      // cas LSL
      else {
	cstart [LSL] = new SMT_CC_Circle (*start_left_forward);
	cend [LSL] = new SMT_CC_Circle (*end_left_backward);
	qi1 [LSL] = new SMT_Configuration (*qa);
	qi2 [LSL] = new SMT_Configuration (*qb);
	length [LSL] = dubins_turn_length (start_left_forward, qa)
	  + configuration_distance (*qa, *qb)
	  + dubins_turn_length (end_left_backward, qb); }
    }

  // cas LSR et sous-cas LS, SR
  if ( internal_mu_tangent_exists (*start_left_forward, *end_right_backward) )
    {
      internal_mu_tangent (*start_left_forward, *end_right_backward, 
			   &qa, &qb);
      // sous-cas LS
      if ( configuration_aligned (*qb, end) ) {
	cstart [LS] = new SMT_CC_Circle (*start_left_forward);
	qi1 [LS] = new SMT_Configuration (*qa);
	length [LS] = dubins_turn_length (start_left_forward, qa)
	  + configuration_distance (*qa, end); }
      // sous-cas SR
      else if ( configuration_aligned (start, *qa) ) {
	cend [SR] = new SMT_CC_Circle (*end_right_backward);
	qi1 [SR] = new SMT_Configuration (*qb);
	length [SR] = configuration_distance (start, *qb)
	  + dubins_turn_length (end_right_backward, qb); }
      // cas LSR
      else {
	cstart [LSR] = new SMT_CC_Circle (*start_left_forward);
	cend [LSR] = new SMT_CC_Circle (*end_right_backward);
	qi1 [LSR] = new SMT_Configuration (*qa);
	qi2 [LSR] = new SMT_Configuration (*qb);
	length [LSR] = dubins_turn_length (start_left_forward, qa)
	  + configuration_distance (*qa, *qb)
	  + dubins_turn_length (end_right_backward, qb); }
    }

  // cas RSL et sous-cas RS, SL
  if ( internal_mu_tangent_exists (*start_right_forward, *end_left_backward) )
    {
      internal_mu_tangent (*start_right_forward, *end_left_backward, 
			   &qa, &qb);
      // sous-cas RS
      if ( configuration_aligned (*qb, end) ) {
	cstart [RS] = new SMT_CC_Circle (*start_right_forward);
	qi1 [RS] = new SMT_Configuration (*qa);
	length [RS] = dubins_turn_length (start_right_forward, qa)
	  + configuration_distance (*qa, end); }
      // sous-cas SL
      else if ( configuration_aligned (start, *qa) ) {
	cend [SL] = new SMT_CC_Circle (*end_left_backward);
	qi1 [SL] = new SMT_Configuration (*qb);
	length [SL] = configuration_distance (start, *qb)
	  + dubins_turn_length (end_left_backward, qb); }
      // cas RSL
      else {
	cstart [RSL] = new SMT_CC_Circle (*start_right_forward);
	cend [RSL] = new SMT_CC_Circle (*end_left_backward);
	qi1 [RSL] = new SMT_Configuration (*qa);
	qi2 [RSL] = new SMT_Configuration (*qb);
	length [RSL] = dubins_turn_length (start_right_forward, qa)
	  + configuration_distance (*qa, *qb)
	  + dubins_turn_length (end_left_backward, qb); }
    }

  // cas RSR et sous-cas RS, SR
  if ( external_mu_tangent_exists (*start_right_forward, *end_right_backward) )
    {
      external_mu_tangent (*start_right_forward, *end_right_backward, 
			   &qa, &qb);
      // sous-cas RS
      if ( configuration_aligned (*qb, end) ) {
	cstart [RS] = new SMT_CC_Circle (*start_right_forward);
	qi1 [RS] = new SMT_Configuration (*qa);
	length [RS] = dubins_turn_length (start_right_forward, qa)
	  + configuration_distance (*qa, end); }
      // sous-cas SR
      else if ( configuration_aligned (start, *qa) ) {
	cend [SR] = new SMT_CC_Circle (*end_right_backward);
	qi1 [SR] = new SMT_Configuration (*qb);
	length [SR] = configuration_distance (start, *qb)
	  + dubins_turn_length (end_right_backward, qb); }
      // cas RSR
      else {
	cstart [RSR] = new SMT_CC_Circle (*start_right_forward);
	cend [RSR] = new SMT_CC_Circle (*end_right_backward);
	qi1 [RSR] = new SMT_Configuration (*qa);
	qi2 [RSR] = new SMT_Configuration (*qb);
	length [RSR] = dubins_turn_length (start_right_forward, qa)
	  + configuration_distance (*qa, *qb)
	  + dubins_turn_length (end_right_backward, qb); }
    }

  // cas LRL
  if ( tangent_circle_exists (*start_left_forward, *end_left_backward) )
    {
      SMT_CC_Circle *middle_right_forward;
      tangent_circle (*start_left_forward, *end_left_backward,
		      &qi1 [LR1L], &qi2 [LR1L], &qi1 [LR2L], &qi2 [LR2L]);

      cstart [LR1L] = new SMT_CC_Circle (*start_left_forward);
      cend [LR1L] = new SMT_CC_Circle (*end_left_backward);
      middle_right_forward = new SMT_CC_Circle (*qi1 [LR1L], 
						false, true, kappa, sigma);
      ci1 [LR1L] = middle_right_forward;
      length [LR1L] = dubins_turn_length (start_left_forward, qi1 [LR1L])
	+ dubins_turn_length (middle_right_forward, qi2 [LR1L])
	+ dubins_turn_length (end_left_backward, qi2 [LR1L]);

      cstart [LR2L] = new SMT_CC_Circle (*start_left_forward);
      cend [LR2L] = new SMT_CC_Circle (*end_left_backward);
      middle_right_forward = new SMT_CC_Circle (*qi1 [LR2L], 
						false, true, kappa, sigma);
      ci1 [LR2L] = middle_right_forward;
      length [LR2L] = dubins_turn_length (start_left_forward, qi1 [LR2L])
	+ dubins_turn_length (middle_right_forward, qi2 [LR2L])
	+ dubins_turn_length (end_left_backward, qi2 [LR2L]);
    }

  // cas RLR
  if ( tangent_circle_exists (*start_right_forward, *end_right_backward) )
    {
      SMT_CC_Circle *middle_left_forward;
      tangent_circle (*start_right_forward, *end_right_backward,
		      &qi1 [RL1R], &qi2 [RL1R], &qi1 [RL2R], &qi2 [RL2R]);

      cstart [RL1R] = new SMT_CC_Circle (*start_right_forward);
      cend [RL1R] = new SMT_CC_Circle (*end_right_backward);
      middle_left_forward = new SMT_CC_Circle (*qi1 [RL1R], 
					       true, true, kappa, sigma);
      ci1 [RL1R] = middle_left_forward;
      length [RL1R] = dubins_turn_length (start_right_forward, qi1 [RL1R])
	+ dubins_turn_length (middle_left_forward, qi2 [RL1R])
	+ dubins_turn_length (end_right_backward, qi2 [RL1R]);

      cstart [RL2R] = new SMT_CC_Circle (*start_right_forward);
      cend [RL2R] = new SMT_CC_Circle (*end_right_backward);
      middle_left_forward = new SMT_CC_Circle (*qi1 [RL2R], 
					       true, true, kappa, sigma);
      ci1 [RL2R] = middle_left_forward;
      length [RL2R] = dubins_turn_length (start_right_forward, qi1 [RL2R])
	+ dubins_turn_length (middle_left_forward, qi2 [RL2R])
	+ dubins_turn_length (end_right_backward, qi2 [RL2R]);
    }
  
 label_fin:  
  // determination du plus court chemin
  cc_dubins_path_type best_path = 
    (cc_dubins_path_type) array_index_min (length, nb_cc_dubins_paths);

  // mesure du temps ecoule
  double time = look_at_chronometer ();

 // contruction du plus court chemin
  SMT_CC_Dubins_Path *path;
  path = new SMT_CC_Dubins_Path (start, end,
				 best_path, kappa, sigma,
				 qi1 [best_path], qi2 [best_path],
				 cstart [best_path], cend [best_path],
				 ci1 [best_path],
				 length [best_path], time);

  // trace des calculs
//    cout << endl << "CC_Dubins" << endl;
//    for (int i = 0; i < nb_cc_dubins_paths; i++) {
//      cout << i << ": "; 
//      if ( !null (qi1[i]) ) { qi1[i]->print (false); }
//      cout << ", "; 
//      if ( !null (qi2[i]) ) { qi2[i]->print (false); }
//      cout << ", " << length[i] << endl;
//    }  
//    cout << "plus court chemin: " << (int) best_path << endl;

  // menage
  delete start_left_forward;
  delete start_right_forward;
  delete end_left_backward;
  delete end_right_backward;
  for (int i = 0; i < nb_cc_dubins_paths; i++) {
    if ( i != best_path ) {
      delete qi1 [i]; delete qi2 [i];
      delete cstart [i]; delete ci1 [i]; delete cend [i];
    }
  }

  //retour du plus court chemin
  return path;
}
