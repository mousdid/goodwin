#ifndef _ODE_SYSTEM_H

#include </home/yassir/Bureau/2A (copie)/PG201/pg201-students-m2-tp3b-mousdid001/libraries/eigen/Eigen/Dense>
#include <fstream>
#include <cmath>
#include <random>

class OdeSystem
{

  protected:
    Eigen::VectorXd _f; // Votre fonction f
    std::ofstream _file_out; // Écriture du fichier
  public:
    // Constructeur par défaut
    OdeSystem();
    // Destructeur par défaut
    virtual ~OdeSystem();
    // Initialiser le nom du fichier solution
    void InitializeFileName(const std::string file_name);
    // Sauvegarde la solution
    virtual void SaveSolution(const double t, const Eigen::VectorXd & sol);
    // Pour récupérer _f
    const Eigen::VectorXd &  GetF()  const;
    // Pour construire _f en fonction de votre système
    virtual void BuildF(const double t, const Eigen::VectorXd & sol)=0;
};


// Classe fille publique d'OdeSystem
class FirstExampleOdeSystem : public OdeSystem
{
public:
void BuildF(const double t, const Eigen::VectorXd & sol); 
//f(X,t) = X

};
class SecondExampleOdeSystem : public OdeSystem
{
public:
void BuildF(const double t, const Eigen::VectorXd & sol); 


};
class ThirdExampleOdeSystem : public OdeSystem
{
public:
void BuildF(const double t, const Eigen::VectorXd & sol); 

};


class LotkaVolterraOdeSystem : public OdeSystem
{
private:
  double _a,_b,_c,_d;  
public:
LotkaVolterraOdeSystem(const double a, const double b, const double c, const double d);//constructeur spécifique

void BuildF(const double t, const Eigen::VectorXd & sol); 


};


class PendulumOdeSystem: public OdeSystem
{
private:
  double _l,_m,_k;  
public:
PendulumOdeSystem(const double m, const double l);//constructeur spécifique 1

PendulumOdeSystem(const double m, const double l,const double k);//constructeur spécifique 2

void BuildF(const double t, const Eigen::VectorXd & sol); 

void SaveSolution(const double t, const Eigen::VectorXd & sol);


};
class DeterministicGoodwinOdeSystem1: public OdeSystem
{
private:
  double  _rho, _beta,  _gamma, _n, _phi, _lambda;  
public:
DeterministicGoodwinOdeSystem1(const double rho, const double beta, const double gamma, const double n, const double phi ,const double lambda);//constructeur spécifique

void BuildF(const double t, const Eigen::VectorXd & sol); 


};

class NonDeterministicOdeSystem : public OdeSystem
{
private:
  double  _rho, _beta,  _gamma, _n, _phi, _lambda; 
std::default_random_engine _generator0;
std::default_random_engine _generator1;
std::normal_distribution<double> _distribution;
public:
NonDeterministicOdeSystem(const double rho, const double beta, const double gamma, const double n, const double phi ,const double lambda);//constructeur spécifique

void BuildF(const double t, const Eigen::VectorXd & sol); 


};




#define _ODE_SYSTEM_H
#endif
