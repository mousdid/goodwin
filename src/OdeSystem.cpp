#ifndef _ODE_SYSTEM_CPP

#include "OdeSystem.h"
#include <iostream>
#include <random>
#include <ctime>



using namespace Eigen;
using namespace std;

// Constructeur par défaut
 OdeSystem::OdeSystem()
{
}

// Destructeur par défaut
OdeSystem::~OdeSystem()
{
}

// Initialisation du nom du fichier
void OdeSystem::InitializeFileName(const std::string file_name)
{
  _file_out.open(file_name);
}

// Renvoie le vecteur _f
const VectorXd & OdeSystem::GetF() const
{
  return this->_f;
}



// Enregistre la solution
// Pour le moment : sol_1, sol_2 ...
void OdeSystem::SaveSolution(const double t, const VectorXd & sol)
{
  this->_file_out << t;
  for (int i = 0 ; i < sol.rows() ; ++i)
  {
    this->_file_out << " " << sol(i);
  }
  this->_file_out << std::endl;
}

// Construit le vecteur f(t, sol)
void FirstExampleOdeSystem::BuildF(const double t, const VectorXd & sol)
{
this->_f = sol; // f(t,X) = X pour résoudre X' = X


}

void SecondExampleOdeSystem::BuildF(const double t, const VectorXd & sol)
{
this->_f.resize(2);
this->_f(0)=-sol(1);
this->_f(1)=sol(0);

}

void ThirdExampleOdeSystem::BuildF(const double t, const VectorXd & sol)
{
this->_f.resize(4);
this->_f=t*sol*sol;
}

LotkaVolterraOdeSystem::LotkaVolterraOdeSystem(const double a, const double b, const double c, const double d)//constructeur spécifique Lotka voltera
{
this->_a=a;
this->_b=b;
this->_c=c;
this->_d=d;
}

void LotkaVolterraOdeSystem::BuildF(const double t, const Eigen::VectorXd & sol)
{
this->_f.resize(2);
this->_f(0)=sol(0)*(this->_a-this->_b*sol(1));
this->_f(1)=sol(1)*(this->_c*sol(0)-this->_d);

}


PendulumOdeSystem::PendulumOdeSystem(const double m, const double l)//constructeur spécifique pendule 1
{
this->_m=m;
this->_l=l;
this->_k=0.;
}

PendulumOdeSystem::PendulumOdeSystem(const double m, const double l, const double k)//constructeur spécifique pendule 2
{
this->_m=m;
this->_l=l;
this->_k=k;
}

void PendulumOdeSystem::BuildF(const double t, const Eigen::VectorXd & sol)
{
const double g=9.81;
this->_f.resize(2);
this->_f(0)=sol(1);
this->_f(1)=(-g/this->_l)*sin(sol(0))-((this->_k)/(this->_m*(pow(this->_l,2))))*sol(1);
}

// Enregistre la solution le cas du pendule
// Pour le moment : sol_1, sol_2 ...
void PendulumOdeSystem::SaveSolution(const double t, const VectorXd & sol)
{
  this->_file_out << t;
  this->_file_out << " " << sol(0);
  this->_file_out << " " << sin(sol(0));
  this->_file_out << " " << -cos(sol(0));
  this->_file_out << std::endl;
}


DeterministicGoodwinOdeSystem1::DeterministicGoodwinOdeSystem1(const double rho, const double beta, const double gamma, const double n, const double phi ,const double lambda)//constructeur spécifique Lotka voltera
{
this->_rho=rho;
this->_beta=beta;
this->_gamma=gamma;
this->_n=n;
this->_phi=phi;
this->_lambda=lambda;
}

void DeterministicGoodwinOdeSystem1::BuildF(const double t, const Eigen::VectorXd & sol)
{
this->_f.resize(2);
this->_f(0)=sol(0)*((1-sol(1))*this->_rho-(this->_beta+this->_gamma+this->_n));
this->_f(1)=sol(1)*((this->_lambda/(1-sol(0)))-(this->_phi+this->_gamma));
}

NonDeterministicOdeSystem::NonDeterministicOdeSystem(const double rho, const double beta, const double gamma, const double n, const double phi ,const double lambda)://constructeur spécifique Lotka voltera
_generator0(123),_generator1 (124),_distribution(0.0,1.0)
{
this->_rho=rho;
this->_beta=beta;
this->_gamma=gamma;
this->_n=n;
this->_phi=phi;
this->_lambda=lambda;


}
void NonDeterministicOdeSystem::BuildF(const double t, const Eigen::VectorXd & sol)
{


this->_f.resize(2);
this->_f(0)=sol(0)*((1-sol(1))*this->_rho-(this->_beta+this->_gamma+this->_n))+0.1*this->_distribution(this->_generator0);
this->_f(1)=sol(1)*(this->_lambda*sol(0)-(this->_phi+this->_gamma))+0.1*this->_distribution(this->_generator1);
}
#define _ODE_SYSTEM_CPP
#endif
