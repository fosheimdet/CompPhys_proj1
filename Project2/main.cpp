
#include <iostream>
#include <armadillo>
#include <iomanip>
#include <string>
#include <chrono>

#include "jacobi.h"
#define CATCH_CONFIG_MAIN //Lets catch define main
#include "catch.hpp"


using namespace arma; 
using namespace std; 
using namespace std::chrono;

vec eigVecCol(4,fill::zeros); // Vector for keeping track of the column indices for the lowest 4 eigenvalues. This is used for the buckling beam problem

double pi = 3.14159;
int n = 100; // Number of  grid points for both the buckling beam- and the harmonic oscillator problem




TEST_CASE("Jacobi vs Armadillo") {
	int N = 10;					//number of gridpoints
	double rho_min = 0; 
	double rho_max = 1; 

	//Initializing all the vectors and matrices needed for the function 'JacobiVsArmadillo'
	mat A = constructA(rho_min, rho_max, N, 0, 0, 0);		
	mat R(N,N, fill::zeros); 
	cx_vec eigVal; 
	cx_mat eigVec;

	double jacobiTime, armaTime; 
	double jacobiSum = 0; 
	double armaSum = 0; 

	//finding the average run time for both algorithms
	for (int i = 0; i < 10; i++) {
		std::tie(jacobiTime, armaTime) = jacobiVsArmadillo(A, R, eigVal, eigVec, eigVecCol, N);  //Breakpoint is set here when timing the algorithms
		jacobiSum += jacobiTime;
		armaSum += armaTime; 
	}

	cout << "Average run time for Jacobi: " << jacobiSum*1e3<< "milli seconds"<< endl; 
	cout << "Average run time for Armadillo: " << armaSum*1e3 << "milli seconds"<<endl; 
}

TEST_CASE("How many similarity transforms are sufficient?") {
	int N = 50; //number of gridpoints
	double h = static_cast<double>(1) / (N + 1);
	double rho_min = 0; 
	double rho_max = 1; 

// Calculating analytical eigenvalues 
	double d = 2 / (h*h);
	double a = -1 / (h*h);
	vec lambda_a(N);
	for (int i = 1; i < N + 1; i++) {
		lambda_a(i - 1) = d + 2 * a*cos((i*pi) / (N + 1));
	}
//Calculating eigenvalues and eigevectors using Jacobi's method
	mat A = constructA(rho_min, rho_max, N, 0, 0, 0); 
	mat R(N, N, fill::zeros);
	jacobiMethod(A, R, N, eigVecCol); 
	int k, l; 
	cout <<"The maximum off-diagonal value is " << maxOffDiag(A, k, l, N) << endl; 
}

TEST_CASE("Harmonic oscillator") {
	double rho_min = 0;
	double rho_max = 5;
	//Performing Jacobi method for harmonic oscillator with no potential
	mat A = constructA(rho_min, rho_max, n, 1, 0, 0.01);
	mat R(n, n, fill::zeros);
	jacobiMethod(A, R, n, eigVecCol);
	writeToFile(R, rho_min, rho_max, n, static_cast<int>(eigVecCol(0)), "non_interacting.txt");

	// Defining arrays containing the different potential strengths and corresponing filenames to store the result
	double omega_array[4] = { 0.01,0.5,1,5 };
	string omega[4] = { "omega0.01.txt", "omega0.5.txt", "omega1.txt", "omega5.txt" };

	//Performing Jacobi's method for the different potential strengths and writing to files
	for (int i = 0; i < 4; i++) {
		mat A = constructA(rho_min, rho_max, n, 1, 1, omega_array[i]);
		mat R(n, n, fill::zeros);
		jacobiMethod(A, R, n, eigVecCol);
		writeToFile(R, rho_min, rho_max, n, static_cast<int>(eigVecCol(0)), omega[i]);
	}
/*	eigVecCol.print();*/ 
}

/////////////UNIT TESTS/////////////////

// Buckling beam. Eigenvectors for the lowest 4 eigenvalues are written to file. 
// The eigenvalues found by Jacobi's method are also tested against those found through Armadillo. 
TEST_CASE("Jacobi's method vs Armadillo"){

	double rho_min = 0; 
	double rho_max = 1; 
	//Setting tolerance for difference between the eigenvalues 
	double eps = 1e-6; 

	// Finding eigenvalues for buckling beam problem through Jacobi's method
	mat A = constructA(rho_min, rho_max, n,0, 0, 0); 
	mat R(n, n, fill::zeros); 
	jacobiMethod(A, R, n, eigVecCol); 
	vec jacobi_eigval = sortEigenvalues(A, n);
	 

	//Using Armadillo to find the eigenvalues 
	cx_vec arma_eig;
	eig_gen(arma_eig, A);
	//Sorting the eigenvalues 
	cx_vec arma_eigval = sort(arma_eig); 

	string eigenmodes[4] = { "bucklingbeam0.txt", "bucklingbeam1.txt", "bucklingbeam2.txt", "bucklingbeam3.txt" };

	//Performing Jacobi's method for the different potential strengths and writing to files
	for (int i = 0; i < 4; i++) {
		mat A = constructA(rho_min, rho_max, n, 0, 0, 0);
		mat R(n, n, fill::zeros);
		jacobiMethod(A, R, n, eigVecCol);
		writeToFile(R, rho_min, rho_max, n, static_cast<int>(eigVecCol(i)), eigenmodes[i]);
	}
	//Comparing the two methods for finding eigenvalues
	for (int i = 0; i < n; i++) {
		REQUIRE( abs(jacobi_eigval(i)-arma_eigval(i)) < eps);
	}
	
}

// Testing that the maxOffDiag function works
TEST_CASE("Maximum non-diagonal element") {
	mat A;
	mat R; R.eye(5, 5); // Setting R to the identity matrix

	//Setting up arbitrary 5x5 symmetric matrix
	A << 5 <<-1 << 0 <<-2 << 0 << endr
	  <<-1 << 4 <<-1 << 0 << 0 << endr
	  << 0 <<-1 << 6 << 0 <<-3 << endr
	  <<-2 << 0 << 0 << 4 << 0 << endr
	  << 0 << 0 <<-3 << 0 << 5 << endr;

	//Comparing the maximum element found through the function "maxOffDiag" vs maximum element found by using armadillo functions. 
	//The matrix is then rotated, and the new max values are compared. 
	int k, l;

	for (int i = 0; i < 10; i++) {
		mat B = abs(A);
		for (int j = 0; j < 5; j++) { //setting diagonal elements to zero in order to exclude then when using A.max()
			B(j, j) = 0; 
	   }
		double maxJacobi = maxOffDiag(A, k, l, 5);
		double maxArma = B.max();
		REQUIRE(maxJacobi == maxArma);
		rotate(A, R, k, l, 5); 
	}

}

TEST_CASE("Eigenvalue test") {
	
	//cout << "Calculating known eigenvalues: " << endl;
	mat A; 
	mat R; R.eye(5, 5);
	vec lambda_a(5);
	//Setting up arbitrary 5x5 symmetric matrix
	A << 5 <<-1 << 0 <<-2 << 0 << endr
	  <<-1 << 4 <<-1 << 0 << 0 << endr
	  << 0 <<-1 << 6 << 0 <<-3 << endr
	  <<-2 << 0 << 0 << 4 << 0 << endr
	  << 0 << 0 <<-3 << 0 << 5 << endr; 
	//Analytical eigenvalues 
	lambda_a(0) = 1.9999999999999982;
	lambda_a(1) = 2.448204782068449;
	lambda_a(2) = 4.099302802201341;
	lambda_a(3) = 6.774340500396769;
	lambda_a(4) = 8.67815191533344;


	double eps = 1.0e-6;

	jacobiMethod(A, R, 5, eigVecCol);
	vec eig = sortEigenvalues(A, 5);

	for (int i = 0; i < 5; i++) {
	/*	cout << (fabs(eig(i) - lambda_a(i))) << endl;*/
		REQUIRE(fabs(eig(i) - lambda_a(i)) < eps);
	}
}



