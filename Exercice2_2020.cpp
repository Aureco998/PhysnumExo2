#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* definir a fonction template pour calculer le produit interne
   entre deux valarray
   inputs:
     array1: (valarray<T>)(N) vecteur de taille N
     array2: (valarray<T>)(N) vecteur de taille N
   outputs:
     produitInterne: (T) produit entre les deux vecteurs
*/
template<typename T> T scalarProduct(valarray<T> const& array1,\
valarray<T> const& array2){
  // compute and return the norm2 of a valarray
  return (array1*array2).sum();
} 

/* definir a fonction template pour calculer la norm2 d'un valarray
   inputs:
     array: (valarray<T>)(N) vecteur de taille N
   outputs:
     norm2: (T) norm2 du vecteur
*/
template<typename T> T norm2(valarray<T> const& array){
  // compute and return the norm2 of a valarray
  return sqrt((array*array).sum());
} 

/* definir a fonction template pour calculer le produit vecteur
   entre 2 valarray de dimension 3
   inputs:
     array1, array2: (valarray<T>)(N) vecteurs de taille N
   outputs:
     produitVecteur: (T) produit vectoriel array1 x aray2 
*/
template<typename T> valarray<T> vectorProduct(valarray<T> const& array1, \
valarray<T> const& array2){
  // initialiser le nouveau valarray
  valarray<T> array3=valarray<T>(3);
  // calculer le produit vecteur
  array3[0] = array1[1]*array2[2] - array1[2]*array2[1]; // premier composante
  array3[1] = array1[2]*array2[0] - array1[0]*array2[2]; // deuxieme composante
  array3[2] = array1[0]*array2[1] - array1[1]*array2[0]; // troisieme composante
  // retourner le array3
  return array3;
} 


/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{

private:
  // definition des constantes
  const double pi=3.1415926535897932384626433832795028841971e0;
  // definition des variables
  double tfin=0.e0;      // Temps final
  unsigned int nsteps=1; // Nombre de pas de temps
  double l_k=0.e0;       // longueur caracteristique du champ magnetique 

  valarray<double> x0=valarray<double>(0.e0,2); // vecteur contenant la position initiale du ballon en
			 		        // utilisant la sequence index-valeur: 0-x, 1-z
  valarray<double> v0=valarray<double>(0.e0,2); // vecteur contenant la vitesse initiale du ballon en
			  		        // utilisant la sequence index-valeur: 0-vx, 1-vz
  unsigned int sampling=1; // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  // TODO calculer l'energie mecanique et le moment magnetique
  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
       write: (bool) ecriture de tous les sampling si faux
  */  
  void printOut(bool write)
  {
    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      double mechanicalEnergy = 1.0/2.0 *mass*(v[0]*v[0] + v[1]*v[1]) + charge* E0 *(x0[1]-x[1]); // TODO calculer l'energie mecanique 
      double magneticMoment   = (mass * (v[0]*v[0] + v[1]*v[1]) / (2*B0)); // TODO calculer le moment magnetique 
      *outputFile << t << " " << x[0] << " " << x[1] << " " \
      << v[0] << " " << v[1] << " " << mechanicalEnergy << " " \
      << magneticMoment << endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // Iteration temporelle, a definir au niveau des classes filles
  virtual void step()=0;

protected:
	//Mis en protected car on en a besoin pour BorisBuneman
	double charge=0.e0; 	 // charge de la particule
	double mass=1.e0;      // mass de la particule
	double B0=0.e0;	 // intensite du champ magnetique
	double E0=0.e0;	 // intensite du champ electrique

  // donnes internes
  double t,dt;        // Temps courant pas de temps
  double chargeOmass; // charge / mass
  valarray<double> x=valarray<double>(2); // Position actuelle de la particule 
  valarray<double> v=valarray<double>(2); // Vitesse actuelle de la particule

   /* Cette methode calcule le champ magnetique en y
     outputs:
       By: (double) champ magnetique en y
  */ 
  double By() const
  {
     return B0*(1.e0+l_k*x[0]);
  }

  /* Cette methode calcule le champ electrique en z
     outputs:
       Ez: (double) champ electrique en z
  */ 
  valarray<double> Ez() const
  {
    valarray<double> E=valarray<double>({0.e0,E0});
    return E;
  }

  // TODO
  /* Cette methode calcule l'acceleration
     outputs:
       a: (valarray<double>)(2) acceleration dans les directions (x,z)
  */
  valarray<double> acceleration(valarray<double>& a) const
  {
    // compute the acceleration
    a[0] = (-(charge/mass)*B0*v[1]) ; // TODO
    a[1] = (charge/mass)*(E0 + B0*v[0]); // TODO
    return a;
  }

public:

  /* Constructeur de la classe Engine
     inputs:
       configFile: (ConfigFile) handler du fichier d'input
  */
  Engine(ConfigFile configFile)
  {
    // variable locale
    
    // Stockage des parametres de simulation dans les attributs de la classe
    tfin     = configFile.get<double>("tfin");		 // lire la temps totale de simulation
    nsteps   = configFile.get<unsigned int>("nsteps");   // lire la nombre de pas de temps
    mass     = configFile.get<double>("mass");		 // lire la mass de la particule
    charge   = configFile.get<double>("charge");	 // lire la charge de la particule
    B0       = configFile.get<double>("B0");		 // lire l intensite du champ magnetique
    l_k      = configFile.get<double>("l_k");		 // lire la longueur du gradient de B
    E0       = configFile.get<double>("E0");		 // lire l intensite du champ electrique
    x0[0]    = configFile.get<double>("x0");		 // lire composante x position initiale
    x0[1]    = configFile.get<double>("z0");		 // lire composante z position initiale
    v0[0]    = configFile.get<double>("vx0");		 // lire composante x vitesse initiale
    v0[1]    = configFile.get<double>("vz0");		 // lire composante z vitesse initiale
    sampling = configFile.get<unsigned int>("sampling"); // lire le parametre de sampling

    dt = tfin / nsteps;          // calculer le time step
    chargeOmass = charge / mass; // coefficient charge/mass

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str()); 
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur virtuel
  virtual ~Engine()
  {
    outputFile->close();
    delete outputFile;
  };

  // Simulation complete
  void run()
  {
    t = 0.e0; // initialiser le temps
    x = x0;   // initialiser la position
    v = v0;   // initialiser la vitesse
    last = 0; // initialise le parametre d'ecriture
    printOut(true); // ecrire premier pas de temps
    for(unsigned int i(0); i<nsteps; ++i) // boucle sur tout pas de temps
    {
      step();  // faire la mise a jour de la simulation 
      printOut(false); // ecrire pas de temps actuel
    }
    printOut(true); // ecrire dernier pas de temps
  };

};

// Extension de la class Engine implementant l'integrateur d'Euler
class EngineEulerImplicite: public Engine
{
private:
  unsigned int maxit=1000; // nombre maximale d iterations
  double tol=1.e12;        // tolerance methode iterative
public:
  
  // construire la class Engine
  EngineEulerImplicite(ConfigFile configFile): Engine(configFile){

    tol = configFile.get<double>("tol"); // lire la tolerance pour la method iterative
    maxit = configFile.get<unsigned int>("maxit"); // lire le nombre d iterations maximale 

  }

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le schema: Euler-Implicit en utilisant la methode de point fixe
     pour resoude le problem implicite
  */
  void step()
  {
    unsigned int iteration=0;
    double error=999e0;
    valarray<double> a=valarray<double>(0.e0,2);
    //valarray<double> xold=valarray<double>(x);
    //valarray<double> vold=valarray<double>(v);
    
    
   
    //Définition de y_n
	valarray<double> y_n({x[0],x[1],v[0],v[1]});
	
   //Définition de y_n+1
	valarray<double> y_n1(y_n);
	
		do{
			
			//Mise à jour de l'accélération
			a = acceleration(a);
			
			//Définition de f(y_n+1)	
			valarray<double> f_y_n1({y_n1[2],y_n1[3],a[0],a[1]});
			
			//Mise à jour des variables
			y_n1 = (y_n + f_y_n1*dt);
			
			 //Définition d^'un vecteur pour l'erreur
			valarray<double> vecterror(0.0, 4);
			vecterror = (y_n1 - y_n - f_y_n1*dt);
			
			error = norm2(vecterror); //mise à jour de l'erreur 
			iteration += 1; //Implementation de l'iteration
			
			//~ //Stockage des valeurs pour la prochaine itération
			//~ xold = x; 
			//~ vold = v;
			
		}while((iteration < maxit) and (error > tol));
		
		x = valarray<double>({y_n1[0], y_n1[1]});
		v = valarray<double>({y_n1[2], y_n1[3]});
		
		t += dt; // mis a jour du temps
  }
};

// Extension de la class Engine implementant l'integrateur d'Euler -> Explicite 
class EngineEuler: public Engine
{
public:
  
  // construire la class Engine
  EngineEuler(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le schema: Euler
  */
  void step()
  {
    valarray<double> a=valarray<double>(0.e0,2);
    x += v*dt; // TODO 
    v += acceleration(a) *dt; // TODO 
    t += dt; // mis a jour du temps
  }
};

// Extension de la class Engine implementant l'integrateur d'Euler-Cromer
class EngineEulerCromer: public Engine
{
public:

  // construire la class Engine
  EngineEulerCromer(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le scheme: Euler-Cromer
  */
  void step()
  {
     valarray<double> a=valarray<double>(0.e0,2);
     
	//Mise à jour de l'accélération
	a = acceleration(a);

     x += v*dt; // TODO
     v[0] += a[0]*dt; // TODO
     v[1] += (charge/mass)*(E0 + B0*v[0])*dt; 

     t += dt; // mis a jour du temps
  }
};

// Extension de la class Engine implementant l'integrateur Runge-Kutta 2
class EngineRungeKutta2: public Engine
{
public:

  // construire la class Engine
  EngineRungeKutta2(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le scheme: Runge-kutta 2
  */
  void step()
  {
    valarray<double> a=valarray<double>(0.e0,2);
    
    // On conserve les valeurs de x et y pour la dernière formule
    valarray<double> xold(x);
	valarray<double> vold(v);

    acceleration(a); //  Mise à jour de l'acceleration
    
    //Définition de k1
	valarray<double> xk1(dt*v);
	valarray<double> vk1(dt*a);
	
	//Mise à jour de x et v
	x += 1.0/2.0 *xk1;
	v += 1.0/2.0 *vk1;
	
	//Mise à jour de l accélération
	acceleration(a);
	
	 //Définition de k2
	valarray<double> xk2(dt*v); // ->(Eq 2.116)
	valarray<double> vk2(dt*a); 
	   
    x = xold + xk2; // TODO
    v = vold + vk2; // TODO

    t += dt; // mis a jour du temps
    
  }
};


// Extension de la class Engine implementant l'integrateur Boris Buneman
class EngineBorisBuneman: public Engine
{

protected:

  // TODO
  // Calcul la rotation des vitesses pour la method de Boris Buneman
  double rotationVitessesBorisBuneman()
  {
	return (charge*B0/ mass);
  }

public:

  // construire la class Engine
  EngineBorisBuneman(ConfigFile configFile): Engine(configFile) {}

  // TODO
  /* Cette methode integre les equations du mouvement en utilisant
     le scheme: Boris-Buneman. 
  */
  void step()
  {
   // Création de la variable beta
	double beta(rotationVitessesBorisBuneman()*dt);
	
	//Définition de x(ti) et v(ti)
	valarray<double> x_(x + (v*dt/2.0));
	valarray<double> v_(v + (charge/mass)*Ez()*dt/2.0);
	
	 //Définition de vx(ti+1) et vz(ti+1)
	double v_xp = (((1-beta*beta/4.0)*v_[0] - beta*v_[1])/(1+beta*beta/4.0)); 
	double v_zp = (((1-beta*beta/4.0)*v_[1] + beta*v_[0])/(1+beta*beta/4.0)); 
	
	//Vecteur v+ 
	valarray<double> vp({v_xp, v_zp});
	
    
    v = (vp + (charge/mass)*Ez()*dt/2.0); // TODO
	x = (x_ + v*dt/2.0); // changement de la position après la vitesse car on a besoin de la vitesse à t+1
	
    t += dt; // mis a jour du temps
    
  }
};

// programme
int main(int argc, char* argv[])
{
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Schema numerique ("Euler"/"E", "EulerCromer"/"EC" ou "RungeKutta2"/"RK2")
  string schema(configFile.get<string>("schema"));

  Engine* engine; // definer la class pour la simulation
  // choisir quel schema utiliser
  if(schema == "Euler" || schema == "E")
  {
    // initialiser une simulation avec schema Euler
    engine = new EngineEuler(configFile);
  }
  else if(schema == "EulerImplicite" || schema == "EI")
  {
    // initialiser une simulation avec schema Euler Implicite
    engine = new EngineEulerImplicite(configFile);
  }
  else if(schema == "EulerCromer" || schema == "EC")
  {
    // initialiser une simulation avec schema Euler-Cromer
    engine = new EngineEulerCromer(configFile);
  }
  else if(schema == "RungeKutta2" || schema == "RK2")
  {
    // initialiser une simulation avec schema runge-kutta 2
    engine = new EngineRungeKutta2(configFile);
  }
  else if(schema == "BorisBuneman" || schema == "BB")
  {
    // initialiser une simulation avec schema Boris Buneman
    engine = new EngineBorisBuneman(configFile);
  }
  else
  {
    cerr << "Schema inconnu" << endl;
    return -1;
  }

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}
