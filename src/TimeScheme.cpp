#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme() : _sys(0)
{
}

// Destructeur
TimeScheme::~TimeScheme()
{
}

// Initialisation de vos différentes variables
void TimeScheme::Initialize(double t0, double dt, VectorXd & sol0, string results, OdeSystem* sys)
{
   this->_dt = dt;
   this->_t = t0 ;
   this->_sol0 = sol0;
   this->_sol = sol0;
   this->_sys = sys;
   if (results.size() > 0)
   {
      this->_sys->InitializeFileName(results);
   }
}

// Schéma en temps par défaut : ici Euler Explicite
// Avancer d'un pas de temps

// Enregistre la solution : fait appel à OdeSystem car la solution
// que l'on souhaite sauvegarder peut être différente de _sol SaveSolution
// le système
void TimeScheme::SaveSolution()
{
   this->_sys->SaveSolution(this->_t, this->_sol);
}

// Renvoie _sol (pratique pour calculer l'ordre de votre méthode)
const VectorXd & TimeScheme::GetIterateSolution() const
{
   return this->_sol;
}


//classes filles

// Schéma en temps : ici Euler Explicite
// Avancer d'un pas de temps

void EulerScheme::Advance()
{ 
   VectorXd ptr_k1;

   this->_sys->BuildF(this->_t, this->_sol);
   ptr_k1=  this->_sys->GetF();
   this->_sol +=this->_dt*ptr_k1;
   this->_t += this->_dt;


}

void RungeKuttaScheme3::Advance()
{
//chaque étape on fait mise a jour de la solution. f contient 
VectorXd ptr_k1, ptr_k2, ptr_k3;

   this->_sys->BuildF(this->_t, this->_sol);
   ptr_k1=  this->_sys->GetF();


   this->_sys->BuildF(this->_t+this->_dt/2., this->_sol+(this->_dt/2.)*(ptr_k1));
   ptr_k2=  this->_sys->GetF();
   

   this->_sys->BuildF(this->_t+this->_dt, this->_sol-this->_dt*(ptr_k1)+2.*this->_dt*(ptr_k2));
   ptr_k3=  this->_sys->GetF();

   this->_sol +=(this->_dt/6.)*ptr_k1+ (this->_dt/6.)*4.*(ptr_k2)+(this->_dt/6.)*(ptr_k3);

  


   this->_t += this->_dt;


}


void RungeKuttaScheme4::Advance()
{
//chaque étape on fait mise a jour de la solution. f contient 
VectorXd ptr_k1, ptr_k2, ptr_k3,ptr_k4;

   this->_sys->BuildF(this->_t, this->_sol);
   ptr_k1=  this->_sys->GetF();


   this->_sys->BuildF(this->_t+this->_dt/2., this->_sol+(this->_dt/2.)*(ptr_k1));
   ptr_k2=  this->_sys->GetF();


   this->_sys->BuildF(this->_t+this->_dt/2., this->_sol+this->_dt/2.*(ptr_k2));
   ptr_k3=  this->_sys->GetF();
  

   this->_sys->BuildF(this->_t+this->_dt, this->_sol+this->_dt*(ptr_k3));
   ptr_k4=  this->_sys->GetF();
   
   this->_sol +=(this->_dt/6.)*ptr_k1+(this->_dt/6.)*2.*(ptr_k2)+ (this->_dt/6.)*2.*(ptr_k3)+ (this->_dt/6.)*(ptr_k4);

  
  


   this->_t += this->_dt;


}




void Explicitmidpoint2::Advance()
{
//chaque étape on fait mise a jour de la solution. f contient 
VectorXd ptr_k;

   this->_sys->BuildF(this->_t, this->_sol);
   ptr_k=  this->_sys->GetF();
   this->_sys->BuildF(this->_t+this->_dt/2., this->_sol+(this->_dt/2.)*ptr_k);
   this->_sol += (this->_dt)*this->_sys->GetF();

   


   this->_t += this->_dt;


}


// AdBashforthSheme3::Initialize_adambashforth(int n, Eigen::VectorXd f_nm1, Eigen::VectorXd f_nm2)
// {

//    this->_n = dt;
//    this->_f_nm1 = f_nm1  ;
//    this->_f_nm2 = f_nm2;
// }
 
void AdBashforthSheme3::Advance()
{
//chaque étape on fait mise a jour de la solution. f contient 
if ((this->_n==1))
{
VectorXd ptr_k1, ptr_k2, ptr_k3;

 this->_sys->BuildF(this->_t, this->_sol);
   ptr_k1=  this->_sys->GetF();


   this->_sys->BuildF(this->_t+this->_dt/2., this->_sol+(this->_dt/2.)*(ptr_k1));
   ptr_k2=  this->_sys->GetF();
   

   this->_sys->BuildF(this->_t+this->_dt, this->_sol-this->_dt*(ptr_k1)+2.*this->_dt*(ptr_k2));
   ptr_k3=  this->_sys->GetF();

   
   this->_f_nm2=ptr_k1;


   this->_sol +=(this->_dt/6.)*ptr_k1+ (this->_dt/6.)*4.*(ptr_k2)+(this->_dt/6.)*(ptr_k3);

  
}
else if ((this->_n==2))
{
VectorXd ptr_k1, ptr_k2, ptr_k3;

 this->_sys->BuildF(this->_t, this->_sol);
   ptr_k1=  this->_sys->GetF();


   this->_sys->BuildF(this->_t+this->_dt/2., this->_sol+(this->_dt/2.)*(ptr_k1));
   ptr_k2=  this->_sys->GetF();
   

   this->_sys->BuildF(this->_t+this->_dt, this->_sol-this->_dt*(ptr_k1)+2.*this->_dt*(ptr_k2));
   ptr_k3=  this->_sys->GetF();

   
   this->_f_nm1=ptr_k1;


   this->_sol +=(this->_dt/6.)*ptr_k1+ (this->_dt/6.)*4.*(ptr_k2)+(this->_dt/6.)*(ptr_k3);

  
}

else
{
   VectorXd f_n;
   this->_sys->BuildF(this->_t, this->_sol);
   f_n=this->_sys->GetF();


   this->_sol += ((this->_dt)/12.)*(23.*this->_sys->GetF()-16.*this->_f_nm1+5.*this->_f_nm2);
   this->_f_nm2= this->_f_nm1;
   this->_f_nm1=f_n;



   


   this->_t += this->_dt;
   this->_n+=1;

}

}


AdBashforthSheme3::AdBashforthSheme3 (int n)
{
   this->_n=n;
}
   


#define _TIME_SCHEME_CPP
#endif
