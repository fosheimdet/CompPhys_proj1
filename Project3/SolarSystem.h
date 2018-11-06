#pragma once

#include "celestialbody.h"
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
using namespace std;

//const double pi = 3.14159;

class SolarSystem {
public:

	vector<CelestialBody> the_bodies;  
	//double the_kineticEnergy; 
	//double the_potentialEnergy; 
	double const GM_sun = 4 * pi*pi;
 
	SolarSystem(); 
	CelestialBody& createCelestialBody(double mass, double xpos, double ypos, double xvel, double yvel, double beta, bool relativistic);
	void calculateForcesAndEnergy(vector<CelestialBody>& the_bodies);
	void IntegrateAndWrite(CelestialBody& theBody, string filename, double t_max, int N, bool use_euler, bool writeToFile);
	int numberOfBodies();
	void printBodies();

	double getKineticEnergy();
	double getPotentialEnergy();
	double getTotalEnergy(); 
	double getAngularMomentum(); 

	void eulerVSVerlet(SolarSystem& theSystem);
	void mercuryPrecession(CelestialBody& theBody, double t_max, int N); 
	




};