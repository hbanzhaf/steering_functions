// *** File: object.h
// *** Author(s): Th. Fraichard
// *** Last modified on 19 Sep 2000 

#ifndef OBJECTS_H
#define OBJECTS_H

// classe generique pour les objets (configurations, chemins)
class SMT_Object 
{
public:  
  // methode pure d'affichage alphanumerique des attributs d'un objet. si
  // eol alors la methode se termine par un retour a la ligne
  virtual void print (bool eol) = 0;
};

#endif
