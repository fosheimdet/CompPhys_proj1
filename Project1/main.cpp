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

int const n = pow(10,2);
//int const n = 10; 


int main() {

	//int n; 
	//cout << "Type the number of grid points you want, n" << endl; 
	//cin >> n; 

	/*cout << problem1e(n,0); */

	//problem1c(n, 1);

	/*mat A = randu<mat>(5, 5);
	vec b = randu<vec>(5);

	A.print("A =");
	b.print("b=");
 
	vec x = solve(A, b);
	x.print("x=");*/
	
	/*arma::vec x = arma::solve(A,b);*/ 


//	Calculates the average CPU time. A breakpoint is used for this. 
	double sum = 0; 
	for (int i = 0; i < 10; i++) {
		sum+=problem1e(n,0);
	}
	double time_avrg = sum / 10; 
	cout << "The average CPU time was: " << time_avrg*pow(10,6) << " microseconds." << endl; 



	return 0; 
}

