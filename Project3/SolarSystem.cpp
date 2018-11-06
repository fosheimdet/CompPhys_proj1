#include "SolarSystem.h"
#include "CelestialBody.h"
#include "integrators.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <math.h>
using namespace std;
using namespace std::chrono;


SolarSystem::SolarSystem() {

}

CelestialBody& SolarSystem::createCelestialBody(double mass, double xpos, double ypos, double xvel, double yvel, double beta, bool relativistic) {
	the_bodies.push_back(CelestialBody(mass,xpos, ypos, xvel, yvel, beta, relativistic));
	return the_bodies.back();    //Return reference to the newly added body
}

int SolarSystem::numberOfBodies(){
	return the_bodies.size();
}

void SolarSystem::calculateForcesAndEnergy(vector<CelestialBody>& the_bodies) {

	for (int i = 0; i < numberOfBodies(); i++) {
		CelestialBody body1 = the_bodies[i];
		for (int j = 0; j < numberOfBodies(); j++) {

			if (j != i) {
				CelestialBody body2 = the_bodies[j];
				double r12 = sqrt(pow((body1.xpos - body2.xpos), 2) + pow((body1.ypos - body2.ypos), 2));
				body1.xacc =body1.xacc - GM_sun * body2.mass*(body1.xpos - body2.xpos) / pow(r12, body1.beta+1);
				body2.xacc =body2.xacc + GM_sun * body2.mass*(body1.xpos - body2.xpos) / pow(r12, body1.beta+1);;

				body1.yacc =body1.yacc - GM_sun * body2.mass*(body1.ypos - body2.ypos) / pow(r12, body1.beta+1);
				body2.yacc = body2.yacc + GM_sun * body2.mass*(body1.ypos - body2.ypos) / pow(r12, body1.beta+1);
				the_bodies[j] = body2;
			}
			
	



		}
		the_bodies[i] = body1; 
/*		the_kineticEnergy += 0.5*body1.mass * (body1.xvel*body1.xvel + body1.yvel*body1.yvel);
		the_potentialEnergy += -GM_sun * body1.mass / body1.r;*/ 
	}
}

void SolarSystem::IntegrateAndWrite(CelestialBody& theBody, string filename, double t_max, int N, bool use_euler, bool writeToFile) {
	ofstream toFile;
	toFile.open(filename, ios::app); 

		calculateForcesAndEnergy(the_bodies);
		if (use_euler) {
			euler(theBody, t_max, N);
		}
		else {
			velocityVerlet(theBody, t_max, N, 1);
		}
		

		if(writeToFile) {
		toFile << setprecision(8) << fixed << theBody.xpos << setw(20) << theBody.ypos << endl;
	}
	


}

void SolarSystem::mercuryPrecession(CelestialBody& theBody, double t_max, int N) {
	ofstream mercuryOrbit;
	mercuryOrbit.open("MercuryRel.txt");
	ofstream precession;
	precession.open("precession.txt");
	int counter = 0; 
	for (int i = 1; i < N; i++) {

		calculateForcesAndEnergy(the_bodies);
		velocityVerlet(theBody, t_max, N, 1);
		mercuryOrbit << setprecision(8) << fixed << theBody.xpos << setw(20) << theBody.ypos << endl; 
		
		if ( theBody.r < 0.390000005) {
			precession << setprecision(8) <<fixed<< counter << setw(20) <<atan(theBody.xpos / theBody.ypos) << endl;
			
			counter++;
		}
	}
 


}


void SolarSystem::printBodies() {
	for (int i = 0; i < the_bodies.size(); i++) {
		cout << "Body " << i << ": " << the_bodies[i].mass << endl; 
	}
}

double SolarSystem::getKineticEnergy() {
	double totalKineticEnergy = 0; 
	for (int i = 0; i < numberOfBodies(); i++) {
	 totalKineticEnergy += 0.5*the_bodies[i].mass * (the_bodies[i].xvel*the_bodies[i].xvel + the_bodies[i].yvel*the_bodies[i].yvel);
	}

	return totalKineticEnergy; 

}

double SolarSystem::getPotentialEnergy() {
	double totalPotentialEnergy = 0; 
	for (int i = 0; i < numberOfBodies(); i++) {
		totalPotentialEnergy += -GM_sun * the_bodies[i].mass / the_bodies[i].r; 
	}
	return totalPotentialEnergy; 
}

double SolarSystem::getTotalEnergy() {
	return getKineticEnergy() + getPotentialEnergy(); 
}

double SolarSystem::getAngularMomentum() {
	double totalAngularMomentum = 0; 

	for (int i = 0; i < numberOfBodies(); i++) {
		totalAngularMomentum += sqrt(pow(the_bodies[i].xpos*the_bodies[i].yvel, 2) + pow(the_bodies[i].ypos*the_bodies[i].xvel, 2)); //Angular momentum of celestial body pr. unit mass
	}
	return totalAngularMomentum; 
}

void SolarSystem::eulerVSVerlet(SolarSystem& theSystem) {

	int tmax = 20; 
	ofstream myFile;
	myFile.open("Euler_vs_Verlet.txt");
	myFile << setprecision(8) << fixed << "N " << setw(20) << "Euler(s) " << setw(20) << "Verlet(s)" << endl;
	
	for (int i = 10; i < 1e5 + 1; i=10 * i) {



		high_resolution_clock::time_point te1 = high_resolution_clock::now();
		theSystem.IntegrateAndWrite(theSystem.the_bodies[0], "Earth.txt", tmax, i, 1,0);
		high_resolution_clock::time_point te2 = high_resolution_clock::now();
		duration<double> eulerTime = duration_cast<duration<double>>(te2 - te1);

		high_resolution_clock::time_point tv1 = high_resolution_clock::now();
		theSystem.IntegrateAndWrite(theSystem.the_bodies[0], "Earth.txt", tmax, i, 0,0);
		high_resolution_clock::time_point tv2 = high_resolution_clock::now();
		duration<double> verletTime = duration_cast<duration<double>>(tv2 - tv1);

		myFile << setprecision(8) << fixed << i << setw(20) << eulerTime.count() << setw(20) << verletTime.count() << endl;
	}
}