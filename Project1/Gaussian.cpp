#include "Gaussian.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std; 

double analyticalSolution(double x);
double ourFunction(double h, int n);
void print(double a[], int l);
void dynamicArrayExample(int n);
void gaussian(int n);


double analyticalSolution(double x) { // Calculates the analytical value of u for a given x.
	return 1 - (1 - exp(-10))*x - exp(-10 * x);
}

double ourFunction(double h, int n) {   // Function to calculculate the elements of f.
	return h * h * 100 * exp(-10 * n*h);
}



void print(double a[], int l) { // Prints arrays. Makes for easy testing of the code. 
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

	double *b = new double[n]; // Allocating memory for all the arrays using dynamic memory.
	double *c = new double[n]; 
	double *a = new double[n]; 

	double *f = new double[n];
	double *v = new double[n]; //Numerical approximation to the analytical solution.
	double *u = new double[n]; //Analytical solution.

	double *x =new double[n+2];  // Initializes the x-array.



	for (int i = 1; i < n + 1; i++) { // Assigns the appropriate values for x given our step size, h. Also assigns values to the vectors b,c and a. 
		x[i] = i * h;
		b[i - 1] = 2;
		c[i - 1] = -1;
		a[i - 1] = -1; 
		/*	cout << x[i] << endl; */
	}


	for (int i = 0; i < n; i++) { // assigns appropriate values for f=h^2f_0(x_i)
		f[i] = ourFunction(h,i + 1);
	}

	for (int i = 1; i < n; i++) { //Performs forward substitution on the matrix A and vector f.
		b[i] = b[i] - (a[i-1]/b[i-1])*c[i-1];
		f[i] = f[i] - (a[i-1]/b[i-1])*f[i-1];
	}


	v[n - 1] = f[n-1]/b[n-1];  // Calculates the last element of v, since v[n] is not included in the array.

	
	for (int i = n - 2; i >= 0; i--) { //Performs backward substitution to find our approximation, v.

		v[i] = (f[i] - c[i]*v[i + 1])/b[i];
	}

	for (int i = 0; i < n; i++) { //Finds the analytical solution, u.
		u[i] = analyticalSolution(x[i + 1]);
	}
	//temporary solution to writing to file. Make a function for it
	ofstream toFile("gaussian_elimination.txt"); 
	

	toFile<<"x" <<setw(20)<<"u"<<setw(20)<<"v" << endl;
	toFile << setprecision(8) << fixed << 0 << setw(20) << 0 << setw(20) << 0 << endl; 
	for (int i = 0; i < n; i++) {
		toFile<<setprecision(8)<<fixed<<x[i+1] <<setw(20)<< u[i]<<setw(20)<<v[i]<< endl; 
	}
	toFile << setprecision(8) << fixed << 1 << setw(20) << 0 << setw(20) << 0 << endl;

	toFile.close(); 

	print(u, n); 
	print(v, n); 

	//Deleting memory allocated to the different vectors
	delete[]b;  
	b = nullptr; 
	delete[]c; 
	c = nullptr; 
	delete[]a; 
	a = nullptr; 

	delete[]f; 
	f = nullptr; 
	delete[]v; 
	v = nullptr; 
	delete[]u; 
	u = nullptr;

	delete[]x; 
	x = nullptr; 
}