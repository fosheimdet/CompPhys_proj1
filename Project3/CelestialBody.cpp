#include "CelestialBody.h"
#include "integrators.h"
#include "SolarSystem.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include <iostream>


using namespace std; 

CelestialBody::CelestialBody(double the_mass, double the_xpos, double the_ypos, double the_xvel, double the_yvel, double the_beta, bool relativistic) { 
	mass = the_mass; 
	xpos = the_xpos; 
	ypos = the_ypos; 
	xvel = the_xvel; 
	yvel = the_yvel;
	beta = the_beta; 

	r = sqrt(xpos*xpos + ypos * ypos); 
	if (!relativistic) {
		xacc = -(4 * pi*pi*xpos) / (pow(r, beta + 1));
		yacc = -(4 * pi*pi*ypos) / (pow(r, beta + 1));
	}
	else {
		xacc = -((4*pi*pi*xpos) / (pow(r, beta + 1)))*(1 + 3 * (pow(xpos*yvel, 2) + pow(ypos*xvel, 2)) / (pow(r, 2)*pow(3e8*60*60*24*365/ 1.496e8, 2)));
		yacc = -((4*pi*pi*ypos) / (pow(r, beta + 1)))*(1 + 3 * (pow(xpos*yvel, 2) + pow(ypos*xvel, 2)) / (pow(r, 2)*pow(3e8 * 60 * 60 * 24 * 365 / 1.496e8, 2)));
	}

}