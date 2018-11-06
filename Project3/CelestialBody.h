#pragma once

#include <string>

const double pi = 3.14159;

using namespace std; 

class CelestialBody {
public: 
	double mass; 
	double r;      //current radius of orbit
	double xpos; 
	double ypos; 
	double xvel; 
	double yvel; 

	double xacc; 
	double yacc; 

	double beta; //Exponent of the radius in the denominator of the force
 

	
	CelestialBody(double the_mass, double the_xpos, double the_ypos, double the_xvel, double the_yvel, double the_beta, bool relativistic);
};