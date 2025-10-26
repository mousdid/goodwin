#ifndef _TIME_SCHEME_H

#include "OdeSystem.h"

class TimeScheme
{
  private:
  // Vecteur initial
   Eigen::VectorXd _sol0;
  protected:
    // Pas de temps
    double _dt;
    // Temps en cours
    double _t;
    // Vecteur initial et vecteur solution
    Eigen::VectorXd  _sol;
    // Pointeur vers le système d'EDO
    OdeSystem* _sys;

  public:
    // Constructeur par défaut
    TimeScheme();
    // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
    // alors le compilateur en génère un implicitement.
    virtual ~TimeScheme();
    // Initialisation de vos différentes variables
    void Initialize(double t0, double dt, Eigen::VectorXd & sol0, std::string name_file, OdeSystem* sys);
    // Enregistre la solution un fichier
    void SaveSolution();
    // Une étape du schéma en temps
    virtual void Advance()=0;
   
    
    // Permet de récupérer _sol
    const Eigen::VectorXd & GetIterateSolution() const;
};

//classes filles

class EulerScheme : public TimeScheme
{
public:
void Advance();

};

class RungeKuttaScheme3 : public TimeScheme
{
public:
void Advance();


};

class RungeKuttaScheme4 : public TimeScheme
{
public:
void Advance();

};
class Explicitmidpoint2 : public TimeScheme
{
public:
void Advance();

};

class AdBashforthSheme3 : public TimeScheme
{
  private:
  int _n;



  public:
  Eigen::VectorXd _f_nm1;
  Eigen::VectorXd _f_nm2;
  void Advance() ;
 
  AdBashforthSheme3 (int n);
};

#define _TIME_SCHEME_H
#endif
