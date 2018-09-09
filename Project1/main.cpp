#include <iostream>
#include <cmath>
#include <fstream>
#include <istream>
#include <string>
#include <vector>
#include <ctime>
#include <chrono>
#include <armadillo>


#include "problems.h"



using namespace std;
using namespace std::chrono;
using namespace arma; 

int const n = pow(10,3);



int main() {

	problem1b(n, 1); 
	//problem1c(n, 1);
	//problem1e(n,0); 

	



//	Calculates the average CPU time. A breakpoint is used for this. 
	//double sum = 0; 
	//for (int i = 0; i < 10; i++) {
	//	sum += problem1b(n, 0);
	//	
	//}
	//double time_avrg = sum / 10; 
	//cout << "The average CPU time was: " << time_avrg*pow(10,6) << " microseconds." << endl; 



	return 0; 
}

