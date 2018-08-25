#include <iostream>
#include <cmath>
using namespace std;

int const n = 4;
double const h = static_cast<double>(1) / (n + 1);

double f(double h,int n); 
void print(double a[],int l);



int main() {

	//double d[8] = { 2, 2, 2, 2, 2, 2, 2, 2 };
	//double e[8] = { -1,-1,-1,-1, -1, -1, -1, -1 };
	//double b[8];
	//double u[8] = { 0 };

	double d[4] = { 2, 2, 2, 2 };
	double e[4] = { -1,-1,-1,-1 };
	double b[4];
	double u[4] = { 0 };

	double x[n+2] = {0};  // Initializes the x-array
	for (int i = 0; i < n + 1; i++) {
		x[i] = 0;
	}

	for (int i = 1; i < n+1 ; i++) { // Assigns the appropriate x values given our step size, h
		x[i] =  i * h; 
	/*	cout << x[i] << endl; */
	}

	for (int i = 0; i < n; i++) { // assigns appropriate values for b=h^2f(x_i)
		b[i] = f(h,i+1); 
	/*	cout << b[i] << endl; */
	}

	for (int i = 1; i < n; i++) {
		d[i] = d[i] - pow(e[i - 1],2) / d[i - 1]; 
		b[i] = b[i] - b[i - 1] * e[i - 1] / d[i - 1]; 
	}

	u[n-1] = b[n-1] / d[n-1]; 

	for (int i = n-2; i >= 0; i--) {

		u[i] = (b[i] - e[i]*u[i + 1]) / d[i]; 
	}
	print(d, n); 
	print(b, n);
	print(u, n); 

	return 0; 
}

double f(double h, int n) {   // Function to calculculate the elements of b
	return h*h*100 * exp(-10 * n*h);
}



void print(double a[], int l) { // Prints arrays
	for (int i = 0; i < l; i++) {
		cout << a[i] << endl; 
	}
	cout << endl; 
}