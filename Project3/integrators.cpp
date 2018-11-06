#include "integrators.h"
#include "CelestialBody.h"
#include "SolarSystem.h"

#include <fstream>
#include <iomanip>

using namespace std;

double GM_sun = 4 * pi*pi;
double vsc0 = (60 * 60 * 24 * 365) / (1.496e8); // Velocity scaling factor from km/s to AU/year


void euler(CelestialBody& theBody, double t_max, int N) {
	//stepsize
	double dt = t_max / N;

	theBody.r = sqrt(theBody.xpos*theBody.xpos + theBody.ypos * theBody.ypos);  //4 FLOPS
	theBody.xacc = -(GM_sun*theBody.xpos) / (pow(theBody.r, 3));  //3FLOPS
	theBody.yacc = -(GM_sun*theBody.ypos) / (pow(theBody.r, 3));  //3FLOPS

	theBody.xpos += theBody.xvel*dt;  //2 FLOPS
	theBody.xvel += theBody.xacc*dt;  //2 FLOPS

	theBody.ypos += theBody.yvel*dt;  //2 FLOPS
	theBody.yvel += theBody.yacc*dt;  //2 FLOPS
 
	//total FLOPPAGE: 18 FLOPS


}

void velocityVerlet(CelestialBody& theBody, double t_max, int N, bool relativistic) {

	    double dt = t_max / N;
	//Finding the position first
		theBody.xpos =theBody.xpos + dt*theBody.xvel + (dt*dt / 2)*theBody.xacc;  //6 FLOPS
		theBody.ypos = theBody.ypos + dt*theBody.yvel + (dt*dt / 2)*theBody.yacc;  //6FLOPS
	
	//Updating r
		theBody.r = sqrt(theBody.xpos*theBody.xpos + theBody.ypos*theBody.ypos); //4FLOPS
		
	//Calculating the acceleration in the next step, and using it to find the velocity
		double next_xacc;
		double next_yacc; 
		if (!relativistic) {
			 next_xacc = -(GM_sun*theBody.xpos) / (pow(theBody.r, theBody.beta + 1));  //3FLOPS
			 next_yacc = -(GM_sun*theBody.ypos) / (pow(theBody.r, theBody.beta + 1));  //3FLOPS
		}
		else {
			 next_xacc = -((GM_sun*theBody.xpos) / (pow(theBody.r, theBody.beta + 1)))*(1+3*(pow(theBody.xpos*theBody.yvel, 2) + pow(theBody.ypos*theBody.xvel, 2))/(pow(theBody.r,2)*pow(3e8*vsc0,2)));
			 next_yacc = -((GM_sun*theBody.ypos) / (pow(theBody.r, theBody.beta + 1)))*(1+3*(pow(theBody.xpos*theBody.yvel, 2) + pow(theBody.ypos*theBody.xvel, 2))/(pow(theBody.r,2)*pow(3e8*vsc0,2)));
		}
	

		theBody.xvel = theBody.xvel + (dt / 2)*theBody.xacc + (dt / 2)*next_xacc; //6FLOPS



		theBody.yvel = theBody.yvel + (dt / 2)*theBody.yacc + (dt / 2)*next_yacc; //6FLOPS

	//Updating acceleration
		theBody.xacc = next_xacc; 
		theBody.yacc = next_yacc; 

		//total FLOPPAGE: 34 FLOPS
}