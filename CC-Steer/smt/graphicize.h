// *** File: graphicize.h
// *** Author(s): Th. Fraichard
// *** Last modified on 9 Jul 2008 

// Description: definition de l'ensemble des fonctions de transformation
// d'objets (heritiers de SMT_Object) en un ensemble d'objets graphiques

#ifndef GRAPHICIZE_H
#define GRAPHICIZE_H

#include <iostream>
using namespace std;

#include <utilities.h>
#include <configuration.h>
#include <paths.h>
#include <display.h>

// ###########################################################################
// fonctions associees
// ###########################################################################

// gestion de la taille par defaut de la demi-branche de la croix
// representant graphiquement une configuration (repere "monde")
void set_cross_size (double size);

// gestion du pas de discretisation par defaut pour l'affichage graphique
// d'une clotoide (repere "monde")
void set_clothoid_step (double step);

// transformation du cc-turn defini par un cc-circle c et une configuration
// d'arrivee q (supposee appartenir au cc-circle et avoir une orientation
// correcte) en un ensemble d'objets graphiques.  le booleen cusps indique
// si l'on se place dans le cas avec ou sans point de rebroussement
void cc_turn_graphicize (SMT_CC_Circle c, SMT_Configuration q, bool cusps,
			 FLTK_Display *display,
			 Fl_Color color, int style, int width);

// ###########################################################################
// fonctions de transformation en un ensemble d'objets graphiques
// ###########################################################################

// une configuration est represente par une croix dont la demi-branche est
// de longueur Cross_Size (unite du "monde").
void graphicize (SMT_Configuration q,
		 FLTK_Display *display,  
		 Fl_Color color = FL_BLACK, 
		 int style = FL_SOLID, int width = 0);

// cc-circle
void graphicize (SMT_CC_Circle c,
		 FLTK_Display *display,  
		 Fl_Color color = FL_BLACK, 
		 int style = FL_SOLID, int width = 0);

// cc-dubins path
void graphicize (SMT_CC_Dubins_Path p,
		 FLTK_Display *display,  
		 Fl_Color color = FL_BLACK, 
		 int style = FL_SOLID, int width = 0);

// cc-reeds-shepp path
void graphicize (SMT_CC_RS_Path p,
		 FLTK_Display *display,  
		 Fl_Color color = FL_BLACK, 
		 int style = FL_SOLID, int width = 0);

#endif
