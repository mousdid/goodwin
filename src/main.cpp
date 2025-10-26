#include "TimeScheme.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
using namespace std;
using namespace Eigen;


int main()
{


	//lecture des donnée d'entrée sur le temps
	std::ifstream entry_file("Data/time_entry.txt");
	double t0, tfinal, dt;
	if (entry_file.is_open()){
		entry_file >> t0 ;
		entry_file >> tfinal;
		entry_file >> dt ;

	}
	else{
		std::cerr << "impossible d'ouvrir le fichier:" << "time_entry.txt" << std::endl;
	}
	entry_file.close();
	// std::cout << t0 << endl;
	// double t0(0.), tfinal(4.), dt(0.005); // temps initial, final, pas de temps
	
	int nb_iterations = int(ceil(tfinal/dt)); // Définition du nombre d'itérations
	dt = tfinal / nb_iterations; // Recalcul de dt
	string results; // nom du fichier résultat
	
	VectorXd sol0, exactSol; // Condition initiale et Solution exacte
	cout << "------------------------------------" << endl;
	
	OdeSystem* sys(0);
	//lecture des parametres du systeme
	std::ifstream parameter_file("Data/parameters.txt");
	double rho,  beta,  gamma,  n,  phi , lambda;
	if (parameter_file.is_open()){
		parameter_file >> rho ;
		parameter_file >> beta;
		parameter_file >> gamma ;
		parameter_file >> n ;
		parameter_file >> phi ;
		parameter_file >> lambda ;

	}
	else{
		std::cerr << "impossible d'ouvrir le fichier:" << "parameters.txt" << std::endl;
	}
	parameter_file.close();
		// cout << rho << beta << gamma << n << phi << lambda;
		
			sys = new  DeterministicGoodwinOdeSystem1( rho,  beta,  gamma,  n,  phi , lambda);
			sol0.resize(2);
			sol0(0) = 0.909; sol0(1) =0.75 ;
			// sol0(0) = (phi+gamma)/lambda; sol0(1) =1-(beta+gamma+n)/rho ;
			results="Goodwin_Results";
int userChoiceScheme; // Choix de l'utilisateur
cout << "------------------------------------" << endl;
	cout << "Choississez le schéma : " << endl;
	cout << "1) Euler Explicite"<< endl;
	cout << "2) Runge Kutta 2" << endl;
	cout << "3) Runge Kutta 3" << endl;
	cout << "4)Runge Kutta 4" << endl;
	cout << "5)Adam Bashforth 3" << endl;
	cin >> userChoiceScheme;
	TimeScheme* time_scheme=NULL;
	auto start_time =chrono::high_resolution_clock::now();
	auto end_time=chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count ();
	switch(userChoiceScheme)	
		{
		case 1:
	time_scheme = new EulerScheme; // Objet de TimeScheme
	results += "_Euler_Explicit.txt";
	time_scheme->Initialize(t0, dt, sol0, results, sys); // Initialisation
	time_scheme->SaveSolution(); // Sauvegarde condition initiale
	start_time = chrono::high_resolution_clock::now();

	for (int n = 0; n < nb_iterations; n++)
	{ // Boucle en temps
	time_scheme->Advance();
	time_scheme->SaveSolution();
	}		
	end_time = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count ();
	cout  << "temps de calcul :" << duration << "milisecondes" << endl;


		break;
		case 2:
	time_scheme = new Explicitmidpoint2; // Objet de TimeScheme
	results += "_RK2.txt";
	time_scheme->Initialize(t0, dt, sol0, results, sys); // Initialisation
	time_scheme->SaveSolution(); // Sauvegarde condition initiale
	start_time = chrono::high_resolution_clock::now();

	for (int n = 0; n < nb_iterations; n++)
	{ // Boucle en temps
	time_scheme->Advance();
	time_scheme->SaveSolution();
	}
	end_time = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count ();
	cout  << "temps de calcul :" << duration << "milisecondes" << endl;

	
		break;
		case 3:
	time_scheme = new RungeKuttaScheme3; // Objet de TimeScheme
	results += "_RK3.txt";
	time_scheme->Initialize(t0, dt, sol0, results, sys); // Initialisation
	time_scheme->SaveSolution(); // Sauvegarde condition initiale
	start_time = chrono::high_resolution_clock::now();

	for (int n = 0; n < nb_iterations; n++)
	{ // Boucle en temps
	time_scheme->Advance();
	time_scheme->SaveSolution();
	}
	end_time = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count ();
	cout  << "temps de calcul :" << duration << "milisecondes" << endl;

		break;
		case 4:
	time_scheme = new RungeKuttaScheme4; // Objet de TimeScheme
	results += "_RK4.txt";
	time_scheme->Initialize(t0, dt, sol0, results, sys); // Initialisation
	time_scheme->SaveSolution(); // Sauvegarde condition initiale
	start_time = chrono::high_resolution_clock::now();

	for (int n = 0; n < nb_iterations; n++)
	{ // Boucle en temps
	time_scheme->Advance();
	time_scheme->SaveSolution();
	}
	end_time = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count ();
	cout  << "temps de calcul :" << duration << "milisecondes" << endl;


		break;
		case 5:
	time_scheme = new AdBashforthSheme3(1); // Objet de TimeScheme
	results += "_AB3.txt";
	time_scheme->Initialize(t0, dt, sol0, results, sys); // Initialisation
	time_scheme->SaveSolution(); // Sauvegarde condition initiale
	start_time = chrono::high_resolution_clock::now();

	for (int n = 1; n < nb_iterations+1; n++)
	

	{ // Boucle en temps


	time_scheme->Advance();
	time_scheme->SaveSolution();
	}
	end_time = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count ();
	cout  << "temps de calcul :" << duration << "milisecondes" << endl;

		break;
		default:
			cout << "Ce choix n'est pas possible ! Veuillez recommencer !" << endl;
		exit(0);
		}


		//Etude de la convergence
		double error1;
		int nb_raffinement=8;
		double pente=0.;
		for (int k = 0; k < nb_raffinement; k++)
		{
		

		VectorXd approxSol1 = time_scheme->GetIterateSolution(); // Temps final
		time_scheme->Initialize(t0, dt/(pow(2.,k+1)), sol0, results, sys);
		for (int n = 0; n < nb_iterations*pow(2.,k+1); n++)
		time_scheme->Advance();
		VectorXd approxSol2 = time_scheme->GetIterateSolution(); // Temps final
		double error = ((((approxSol1-approxSol2)).array().abs()).sum());
		if (error<pow(10.,-12))
		{break;}
		cout << "------------------------------------" << endl;
		cout << "Erreur = " << error<< " entre les pas " << dt/pow(2.,k) << " et "<< dt/pow(2.,k+1) << endl;
		if (k>0) 
		{
		pente=log2(error1/error);
		cout << "pente = " << pente << " entre les pas " << dt/pow(2.,k) << " et "<< dt/pow(2.,k+1) << endl;}
		
		
		error1=error;
		}
		cout << "------------------------------------" << endl;
		cout << "Ordre de la méthode = " << pente << endl;
	

delete time_scheme;




































	delete sys;
	
return 0;
}