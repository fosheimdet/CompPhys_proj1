#include "Gaussian.h"
#include <iostream>
#include <cmath>

using namespace std; 

double analyticalSolution(double x);
double f(double h, int n);
void print(double a[], int l);
void dynamicArrayExample(int n);
void gaussian(int n);


double analyticalSolution(double x) { // Calculates the analytical value of u for a given x
	return 1 - (1 - exp(-10))*x - exp(-10 * x);
}

double f(double h, int n) {   // Function to calculculate the elements of b
	return h * h * 100 * exp(-10 * n*h);
}



void print(double a[], int l) { // Prints arrays
	for (int i = 0; i < l; i++) {
		cout << a[i] << endl;
	}
	cout << endl;
}


void dynamicArrayExample(int n) {

	double *numbers = new double[n]; 
	
	for (int i = 0; i < n; i++) {
		cout << "Input number: ";
		cin >> numbers[i];
	}
		cout << "You entered: "; 
		for (int i = 0; i < n; i++) {
			cout << numbers[i] << " "; 
		}
		cout << endl; 

		delete[] numbers; 
		numbers = nullptr; 

	
}


void gaussian(int n) {

	double const h = static_cast<double>(1) / (n + 1);

	double *d = new double[n]; // Allocating memory for all the arrays, using dynamic memory
	double *e = new double[n]; 
	double *b = new double[n];
	double *v = new double[n]; //Numerical approximation to the analytical solution
	double *u = new double[n]; //Analytical solution

	double *x =new double[n+2];  // Initializes the x-array


	for (int i = 1; i < n + 1; i++) { // Assigns the appropriate values for x given our step size, h. Also assigns values to the vectors d and e. 
		x[i] = i * h;
		d[i - 1] = 2;
		e[i - 1] = -1;
		/*	cout << x[i] << endl; */
	}

	for (int i = 0; i < n; i++) { // assigns appropriate values for b=h^2f(x_i)
		b[i] = f(h, i + 1);
		/*	cout << b[i] << endl; */
	}

	for (int i = 1; i < n; i++) { //Performs forward substitution on the matrix A and vector b
		d[i] = d[i] - pow(e[i - 1], 2) / d[i - 1];
		b[i] = b[i] - b[i - 1] * e[i - 1] / d[i - 1];
	}

	v[n - 1] = b[n - 1] / d[n - 1];  // Calculates the last element manually, since u_num[n] is not included in the array u_num

	for (int i = n - 2; i >= 0; i--) { //Performs backward substitution to find our approximation for u

		v[i] = (b[i] - e[i] * v[i + 1]) / d[i];
	}

	for (int i = 0; i < n; i++) { //Finds the analytical solution, u
		u[i] = analyticalSolution(x[i + 1]);
	}

	print(u, n); 
	print(v, n); 

	delete[]d; 
	d = nullptr; 
	delete[]e; 
	e = nullptr; 
	delete[]b; 
	b = nullptr; 
	delete[]v; 
	v = nullptr; 
	delete[]u; 
	u = nullptr; 
	delete[]x; 
	x = nullptr; 
}