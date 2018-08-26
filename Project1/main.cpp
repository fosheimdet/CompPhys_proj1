#include <iostream>
#include <cmath>


using namespace std;

int const n = 16;
double const h = static_cast<double>(1) / (n + 1);

double analyticalSolution(double x);
double f(double h,int n); 
void print(double a[],int l);



int main() {

	//double d[8] = { 2, 2, 2, 2, 2, 2, 2, 2 };
	//double e[8] = { -1,-1,-1,-1, -1, -1, -1, -1 };
	//double b[8];
	//double u[8] = { 0 };

	double d[n] = { 0 }; 
	double e[n] = { 0 }; 
	double b[n] = { 0 }; 
	double u_num[n] = { 0 };
	double u[n] = { 0 }; 
	

	double x[n+2] = {0};  // Initializes the x-array


	for (int i = 1; i < n+1 ; i++) { // Assigns the appropriate values for x given our step size, h. Also assigns values to the vectors d and e. 
		x[i] =  i * h; 
		d[i - 1] = 2; 
		e[i - 1] = -1;
	/*	cout << x[i] << endl; */
	}

	for (int i = 0; i < n; i++) { // assigns appropriate values for b=h^2f(x_i)
		b[i] = f(h,i+1); 
	/*	cout << b[i] << endl; */
	}

	for (int i = 1; i < n; i++) { //Performs forward substitution on the matrix A and vector b
		d[i] = d[i] - pow(e[i - 1],2) / d[i - 1]; 
		b[i] = b[i] - b[i - 1] * e[i - 1] / d[i - 1]; 
	}

	u_num[n-1] = b[n-1] / d[n-1];  // Calculates the last element manually, since u_num[n] is not included in the array u_num

	for (int i = n-2; i >= 0; i--) { //Performs backward substitution to find our approximation for u

		u_num[i] = (b[i] - e[i]*u_num[i + 1]) / d[i]; 
	}

	for (int i = 0; i < n; i++) { //Finds the analytical solution, u
		u[i] = analyticalSolution(x[i+1]); 
	}
	

	//print(d, n); 
	//print(b, n);
	print(u, n);
	print(u_num, n); 

	return 0; 
}

double analyticalSolution(double x) {
	return 1 - (1 - exp(-10))*x - exp(-10 * x);
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