#include <iostream>	
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>


#include "integrators.h"
#include "CelestialBody.h"
#include "SolarSystem.h"

using namespace std; 
using namespace std::chrono;

//Setting number of integration points and maximum time for the solar system
int N_system = 10000; 
double tmax_system = 200;

//Setting number of integration points and maximum time for the relativistic orbit of Mercury
int N_mercury = 1000000; 
double tmax_mercury = 100; 


double vsc = (60 * 60 * 24 * 365) / (1.496e8); // Velocity scaling factor from km/s to au/year

int main() {

	// Using the default constructor to create an instantiation of SolarSystem
	SolarSystem solarsystem; 

	//Adding the planets of the solar system

/*	solarsystem.createCelestialBody(1, -1 / (1 + 3e6), 0, 0, (-(3e-6)*29.79-(9.5e-4)*13.07)*vsc, 2);*/ //The sun
	solarsystem.createCelestialBody(1.6425e-7, 0.39, 0, 0, 48 * vsc, 2, 0); //Mercury
	solarsystem.createCelestialBody(2.4335e-6, 0.72, 0, 0, 35 * vsc,2,0); //Venus
	solarsystem.createCelestialBody(3e-6, 1, 0, 0, 29.79*vsc, 2, 0);  //Earth 
	solarsystem.createCelestialBody(3.195e-7, 1.524, 0, 0, 24.1*vsc,2,0); //Mars
	solarsystem.createCelestialBody(9.5e-4, 5.2, 0, 0, 13.07*vsc,2,0);  //Jupiter
	solarsystem.createCelestialBody(2.75e-4, 9.54, 0, 0, 9.6*vsc, 2,0); //Saturn 
	solarsystem.createCelestialBody(4.4e-5, 19.19, 0, 0, 6.80*vsc, 2,0); //Uranus
	solarsystem.createCelestialBody(5.15e-5, 30.06, 0, 0, 5.43*vsc, 2,0); //Neptune
	
	//removing existing files, since we append to files. 
	remove("Mercury.txt"); 
	remove("Venus.txt"); 
	remove("Earth.txt"); 
	remove("Mars.txt"); 
	remove("Jupiter.txt"); 
	remove("Saturn.txt"); 
	remove("Uranus.txt"); 
	remove("Neptune.txt");



	for (int i = 1; i < N_system; i++) {
		solarsystem.IntegrateAndWrite(solarsystem.the_bodies[0], "Mercury.txt", tmax_system, N_system, 0, 1);
		solarsystem.IntegrateAndWrite(solarsystem.the_bodies[1], "Venus.txt", tmax_system, N_system, 0, 1);
		solarsystem.IntegrateAndWrite(solarsystem.the_bodies[2], "Earth.txt", tmax_system, N_system, 0, 1);
		solarsystem.IntegrateAndWrite(solarsystem.the_bodies[3], "Mars.txt", tmax_system, N_system, 0, 1);
		solarsystem.IntegrateAndWrite(solarsystem.the_bodies[4], "Jupiter.txt", tmax_system, N_system, 0, 1);
		solarsystem.IntegrateAndWrite(solarsystem.the_bodies[5], "Saturn.txt", tmax_system, N_system, 0, 1);
		solarsystem.IntegrateAndWrite(solarsystem.the_bodies[6], "Uranus.txt", tmax_system, N_system, 0, 1);
		solarsystem.IntegrateAndWrite(solarsystem.the_bodies[7], "Neptune.txt", tmax_system, N_system, 0, 1);

  }
	solarsystem.eulerVSVerlet(solarsystem);    //Timing the Euler and Verlet algorithms

	SolarSystem sun_mercury;
	sun_mercury.createCelestialBody(1.6425e-7, 0.39, 0, 0, 48 * vsc, 2, 1);
	sun_mercury.mercuryPrecession(solarsystem.the_bodies[0], tmax_mercury, N_mercury); 


	system("pause"); 
}