
#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>

#include "jacobi.h"

using namespace std; 
using namespace arma; 
using namespace std::chrono; 


arma::mat constructA(double rho_min, double rho_max, int n, bool potential, bool interacting,double omega) {
	double h = static_cast<double>(rho_max-rho_min) / (n + 1);

	mat A(n, n, fill::zeros);

	if (potential != true) {

		vec d(n, fill::zeros); //Initializing diagonal elements
		double e = -1 / (h*h); // Off diagonal elements
		for (int i = 0; i < n; i++) {
			double rho = (i + 1)*h;
			d(i) = 2 / (h*h) ; // Diagonal elements
		}

		for (int i = 0; i < n - 1; i++) {
			A(i, i) = d(i);
			A(i, i + 1) = e;
			A(i + 1, i) = e;
		}
		A(n - 1, n - 1) = d(n - 1);

	}

	if (potential == true && interacting != true) { //Tridiagonal toeplitz matrix for no potential
		vec d(n, fill::zeros); //Initializing diagonal elements
		double e = -1 / (h*h); // Off diagonal elements
		for (int i = 0; i < n; i++) {
			double rho =(i+1)*h;
			d(i) = 2 / (h*h) + rho * rho; // The potential is included here by the rho term
		}

		for (int i = 0; i < n - 1; i++) {
			A(i, i) = d(i);
			A(i, i + 1) = e;
			A(i + 1, i) = e;
		}
		A(n - 1, n - 1) = d(n - 1);


	}
	if(potential == true && interacting == true) { //If potential = true, do the following
		vec d(n, fill::zeros); //Initializing diagonal elements
		double e = -1 / (h*h); // Off diagonal elements
		for (int i = 0; i < n; i++) {
			double rho = (i+1)*h;
			
			d(i) = 2 / (h*h) + omega*omega*rho * rho + 1/(rho); // The potential is included here by the rho term
		}
		for (int i = 0; i < n - 1; i++) {
			A(i, i) = d(i);
			A(i, i + 1) = e;
			A(i + 1, i) = e;
		}
		A(n - 1, n - 1) = d(n - 1);
	}
	return A; 
}



void jacobiMethod(arma::mat &A, arma::mat &R, int n, vec &eigVecCol) {
	// Setting up the eigenvector matrix 

	for (int i = 0; i < n; i++) {
		R(i, i) = 1;
	}

	int k, l;
	double epsilon = 1.0e-8;  //Setting max value for off-diag elements
	double max_iter = n*n*n;   //Maximum number of Jacobi rotations
	int iterations = 0;			// Counts the number of iterations
	double max_offdiag = maxOffDiag(A, k, l, n);  // Finding the max off-diag value

	//Performing the amount of rotations necessary for epsilon > max_offdiag
	while (max_offdiag > epsilon && iterations <= max_iter) {
		rotate(A, R, k, l, n);
		max_offdiag = maxOffDiag(A, k, l, n);
		iterations++;
	}

	//Finding the column indices for the 4 lowest eigenvalues
	/*vec eigVecCol(4);*/  // Vector for keeping track of the column indices for the lowest 4 eigenvalues. This is used for the buckling beam problem
	int ground = 0; 
	int first = 0; 
	int second = 0; 
	int third = 0; 


	for (int i = 0; i < n; i++) {
		if (A(ground, ground) > A(i, i)) {
			ground = i;
		}
		eigVecCol(0) = ground;
	}
	for (int i = 0; i < n; i++) {
		if (A(first, first) > A(i, i) && A(i,i)> A(ground, ground)) {
			first = i;
		}
		eigVecCol(1) = first;
	}
	for (int i = 0; i < n; i++) {
		if (A(second, second) > A(i, i) && A(i,i) > A(first, first)) {
			second = i;
		}
		eigVecCol(2) = second;
	}
	for (int i = 0; i < n; i++) {
		if (A(third, third) > A(i, i)  && A(i,i) > A(second, second)) {
			third = i;
		}
		eigVecCol(3) = third;
	}
}

double maxOffDiag(arma::mat & A, int &k, int &l, int n) {

	double maximumOff = 0; // Maximum value off the diagonal

	for (int i = 0; i < n; i++) {
		//Only include upper triangular part, since A is symmetric. 
		for (int j = i + 1; j < n; j++) {

			double aij = fabs(A(i, j)); 
			if (aij > maximumOff) {
				//Update maximum value
				maximumOff = aij; 
				// Update position
				k = i; 
				l = j; 
			}

		}
	}
	return maximumOff; 
}

void rotate(arma::mat &A, arma::mat &R, int k, int l, int n) {

	//Defining tau
	double t, tau; 
	tau = (A(l, l) - A(k, k)) / (2 * A(k, l)); 

	// Calculating t=s/c as a function of tau. 
	// t can be either positive or negative, depending on the sign 
	// in the expression for t. 

	
	if (tau >= 0) {
		t = 1.0 / (tau + sqrt(1 + tau * tau));   // tan(theta)
	}
	else {
		t = -1.0 / (-tau + sqrt(1 + tau * tau)); // tan(theta)
	}

	double c = 1.0 / (sqrt(1 + t * t)); // cos(theta)
	double s = c * t;					// sin(theta)

	if (A(k, l) == 0) {
		c = 1.0; 
		s = 0.0; 
	}

	// Performing rotation

	double a_kk, a_ll, a_ik, a_il, r_ik, r_il; 
	a_kk = A(k, k); 
	a_ll = A(l, l); 
	//Chaning the matrix elements with indices k and l
	A(k, k) = c * c* a_kk - 2.0*c*s*A(k, l) + s * s*a_ll; 
	A(l, l) = s * s*a_kk + 2 * c*s*A(k, l) + c * c*a_ll; 
	A(k, l) = 0.0; // Hard-coding of zeros
	A(l, k) = 0.0; 

	// Changing the remaining elements 
	for (int i = 0; i < n; i++) {
		if (i != k && i != l) {
			a_ik = A(i, k); 
			a_il = A(i, l); 
			A(i, k) = c * a_ik - s * a_il; 
			A(k, i) = A(i, k);  //Due to symmetry
			A(i, l) = c * a_il + s * a_ik; 
			A(l, i) = A(i, l); 
		}
		//Computing the eigenvectors
		r_ik = R(i, k); 
		r_il = R(i, l); 
		R(i, k) = c * r_ik - s * r_il; 
		R(i, l) = c * r_il + s * r_ik; 
	}

	return; 



}


arma::vec sortEigenvalues(arma::mat A, int n) {
	vec eig(n); 
	//Placing the diagonals in the vector eig
	for (int i = 0; i < n; i++) {
		eig(i) = A(i, i); 
	}
	// Sorting the eigenvalues 
	vec sorted = sort(eig);
	return sorted;



}

void writeToFile(arma::mat &R, double rho_min, double rho_max, int n, int minIndex, string filename) {

	//Setting up rho vector. This could be done in python, but one would not be able to alter rho_min and rho_max. 
	double h = static_cast <double> (rho_max - rho_min) / (n + 1);
	vec rho(n);
	for (int i = 0; i < n; i++) {
		rho(i) = (i + 1)*h;
	}

	//Writing rho and the eigenvector for the groundstate to file
	ofstream toFile(filename);
	toFile << "rho" << setw(20) << "u" << endl;
	toFile << setprecision(8) << fixed << rho_min << setw(20) << 0 << endl; 
	for (int i = 0; i < n; i++) {
		toFile << setprecision(8) << fixed << rho(i) << setw(20) << R(i, minIndex) << endl;
	}
	toFile << setprecision(8) << fixed << rho_max << setw(20) << 0 << endl;

	toFile.close();
}

// The following function computes the eigenvalues and eigenvectors of A using Jacobi's method
// and the Armadillo library. The function returns the run time for both algorithms. 
std::tuple<double,double> jacobiVsArmadillo(arma::mat &A, arma::mat& R, arma::cx_vec eigVal, arma::cx_mat eigVec,arma::vec &eigVecCol, int n) {
	
	high_resolution_clock::time_point tj1 = high_resolution_clock::now();
	jacobiMethod(A, R, n, eigVecCol); 
	high_resolution_clock::time_point tj2 = high_resolution_clock::now();
	duration<double> jacobiTime = duration_cast<duration<double>>(tj2 - tj1);

	high_resolution_clock::time_point ta1 = high_resolution_clock::now();
	eig_gen(eigVal, eigVec, A);
	high_resolution_clock::time_point ta2 = high_resolution_clock::now();
	duration<double> armaTime = duration_cast<duration<double>>(ta2 - ta1);

	return std::make_tuple(jacobiTime.count(), armaTime.count());
}