#include "problems.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <istream>
#include <string>
#include <vector>
#include <ctime>
#include <chrono>
#include <armadillo>

#include <iomanip>


using namespace std; 
using namespace std::chrono;
using namespace arma; 

/*double analyticalSolution(double x);
double ourFunction(double h, int n);
void print(double a[], int l);
void problem1b(int n);
void problem1c(int n);*/ 


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

// Arguments for problem1b and 1c are numbers of gridpoints and a bool for deciding wether to write to file or not.
// The functions return CPU time for forward and backward substitution
double problem1b(int n,bool writeToFile) {

	


	double const h = static_cast<double>(1) / (n + 1);
	//cout << "Step size h: " << h << endl;

	double *a = new double[n];
	double *b = new double[n]; // Allocating memory for all the arrays using dynamic memory.
	double *c = new double[n]; 


	double *f = new double[n];
	double *v = new double[n]; //Numerical approximation to the analytical solution.
	double *u = new double[n]; //Analytical solution.
	double *e = new double[n]; //Error

	double *x =new double[n+2];  // Initializes the x-array.
	


	for (int i = 1; i < n + 1; i++) { // Assigns the appropriate values for x given our step size, h. Also assigns values to the vectors b,c and a. 
		x[i] = i * h;
	}


	//////////////////////////////////////////////////////////////
    // CODE TO READ VECTORS a,b AND c FROM FILE/THROUGH CONSOLE //
	//////////////////////////////////////////////////////////////

	// The file "matrixdata.txt" contains the elements of a, b and c stored in separate columns. 

	const int COLUMNS = 3;

	vector <vector <double>> data; // Creates a vector to contain a, b and c. 
	string filename = "matrixdata232.txt";
	ifstream ifile(filename);

	if (ifile.is_open()) {
		int num;

		vector <double> numbers_in_line; //Creates a vector containing one row of our file 

		string dummyline;
		int skiplines = 4; 
		for (int i = 0; i < skiplines; i++) {
			getline(ifile, dummyline);
		}
		

		for(int i = 0; i<3*n; i++) { //Reads the contents of the file and allocates each row of the file to a vector in "data".
			ifile >> num; 
			numbers_in_line.push_back(num);
			if (numbers_in_line.size() == COLUMNS) {
				data.push_back(numbers_in_line);
				numbers_in_line.clear();
			}
		}
	}
	else {
		
		//cout << "Error opening the file. Assigning values for a, b and c automatically.  \n";

		for (int i = 0; i < n; i++) { // Assigns values for a,b and c if the file was not found.
			a[i] = -1; 
			b[i] = 2; 
			c[i] = -1; 
		}
		
	}

	
	
	
	for (int i = 0; i < data.size(); i++) {//The following code assigns the appropriate elements of data to a, b and c. 
		a[i] = data[i][0];
		b[i] = data[i][1]; 
		c[i] = data[i][2]; 

	}

	ifile.close(); 


	//Lets the user input a, b and c from the console

	//cout << "type the elements of a,b, and c vector" << endl; 
	//for (int i = 0; i < n; i++) {
	//	cout << "Type the " << i << "th element of a" << endl; 
	//	cin >> a[i];
	//	cout << "Type the " << i << "th element of b" << endl;
	//	cin >> b[i];
	//	cout << "Type the " << i << "th element of c" << endl;
	//	cin >> c[i];

	//}

	//////////////////////////////////////////////////////////////

	

	for (int i = 0; i < n; i++) { // assigns appropriate values for f=h^2f_0(x_i)
		f[i] = ourFunction(h,i + 1);
	}
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	//Performs forward substitution on the matrix A and vector f.
	for (int i = 1; i < n; i++) { 
		double ratio = a[i - 1] / b[i - 1];  // 1 FLOP
		b[i] = b[i] - (ratio)*c[i-1];		// 2 FLOPS
		f[i] = f[i] - (ratio)*f[i-1];		// 2 FLOPS
	}


	v[n - 1] = f[n-1]/b[n-1];  // Calculates the last element of v, since v[n] is not included in the array.

	
	for (int i = n - 2; i >= 0; i--) { //Performs backward substitution to find our approximation, v.

		v[i] = (f[i] - c[i]*v[i + 1])/b[i];    //3 FLOPS
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	/*cout << duration_cast<chrono::milliseconds>(t1 - t2).count()<< endl;*/
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	

	for (int i = 0; i < n; i++) { //Finds the analytical solution, u.
		u[i] = analyticalSolution(x[i + 1]);
	}


	//calculate error
	for (int i = 0; i < n; i++) {
		e[i] = log10(abs((v[i] - u[i]) / u[i]));
	}

	//finding maximum value in error array

	double max = e[0];
	for (int i = 0; i < n; i++) {
		if (abs(e[i]) > abs(max)) max = e[i];
	}
	cout << "The maximum relative error was " << max << endl;


	//Writes to file if writeToFile == 1
	if (writeToFile) {
		ofstream toFile("gaussian_elimination.txt");


		toFile << "x" << setw(20) << "u" << setw(20) << "v" << setw(20) << "e" << endl;
		/*error << "e" << endl;*/
		toFile << setprecision(8) << fixed << 0 << setw(20) << 0 << setw(20) << 0 << setw(20) << 0 << endl;
		for (int i = 0; i < n; i++) {
			toFile << setprecision(8) << fixed << x[i + 1] << setw(20) << u[i] << setw(20) << v[i] << setw(20) << e[i] << endl;

		}
		toFile << setprecision(8) << fixed << 1 << setw(20) << 0 << setw(20) << 0 << setw(20) << 0 << endl;

		toFile.close();

	}


	//Deleting memory allocated to the different vectors
	delete[]b;   
	delete[]c; 
	delete[]a; 
 

	delete[]f; 
	delete[]v; 
	delete[]u;
	delete[]e; 
;
	delete[]x;  
	return time_span.count();

}

double problem1c(int n, bool writeToFile) {

	

	double const h = static_cast<double>(1) / (n + 1);
	

	double *b = new double[n]; // Allocating memory for all the arrays using dynamic memory.



	double *f = new double[n];
	double *v = new double[n]; //Numerical approximation to the analytical solution.
	double *u = new double[n]; //Analytical solution.
	double *e = new double[n]; //Error

	double *x = new double[n + 2];  // Initializes the x-array.



	for (int i = 1; i < n + 1; i++) { // Assigns the appropriate values for x given our step size, h. Also assigns values to the vectors b,c and a. 
		x[i] = i * h;
		b[i - 1] = 2;
		//c[i - 1] = -1;
		//a[i - 1] = -1;
		/*	cout << x[i] << endl; */
	}

	

	for (int i = 0; i < n; i++) { // assigns appropriate values for f=h^2f_0(x_i)
		f[i] = ourFunction(h, i + 1);
	}
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for (int i = 1; i < n; i++)
	{ //Performs forward substitution on the matrix A and vector f.
		b[i] = b[i] - 1/b[i - 1];              //2 FLOPS
		f[i] = f[i] + (f[i - 1] / b[i - 1]); //2 FLOPS
	}


	v[n - 1] = f[n - 1] / b[n - 1];  // Calculates the last element of v, since v[n] is not included in the array.


	for (int i = n - 2; i >= 0; i--) { //Performs backward substitution to find our approximation, v.

		v[i] = (f[i] + v[i + 1])/b[i];  // 2 FLOPS
	}

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	//cout << "Symmetric tridiagonal matrix: " << time_span.count() << " seconds";
	//cout << endl;

	for (int i = 0; i < n; i++) { //Finds the analytical solution, u.
		u[i] = analyticalSolution(x[i + 1]);
	}

	//calculate error
	for (int i = 0; i < n; i++) {
		e[i] = log10(abs((v[i] - u[i]) / u[i]));
	}

	//finding maximum value in error array

	double max = e[0];
	for (int i = 0; i < n; i++) {
		if (abs(e[i]) > abs(max)) max = e[i];
	}
	cout <<"The maximum relative error was "<< max << endl; 
	
	
	if (writeToFile) {
		ofstream toFile("gaussian_elimination.txt");


		toFile << "x" << setw(20) << "u" << setw(20) << "v" << endl;
		toFile << setprecision(8) << fixed << 0 << setw(20) << 0 << setw(20) << 0 << endl;
		for (int i = 0; i < n; i++) {
			toFile << setprecision(8) << fixed << x[i + 1] << setw(20) << u[i] << setw(20) << v[i] << endl;
		}
		toFile << setprecision(8) << fixed << 1 << setw(20) << 0 << setw(20) << 0 << endl;

		toFile.close();
	}

	//Deleting memory allocated to the different vectors
	delete[]b;


	delete[]f;
	delete[]v;
	delete[]u;
	delete[]e;  

	delete[]x;

	return time_span.count();


}

double problem1e(int n, bool writeToFile) {

	double h = static_cast<double>(1) / (n + 1);

	mat A(n+2, n+2, fill::zeros);

	vec a(n+2); a.fill(-1);
	vec b(n+2); b.fill(2);
	vec c(n+2); c.fill(-1);
	vec x(n+2);
	vec f(n+2); 
	vec v(n+2);
	/*v[0] = 0; 
	v[n + 2] = 1; */





	for (int i = 0; i < n+1; i++) {
		A(i, i) = b(i);
		A(i, i + 1) = a(i);
		A(i + 1, i) = c(i);
	}
	A(n+1, n+1) = b(n+1);
	//A.print();

	for (int i = 0; i < n+2; i++) {
		x(i) = i * h; 
		f(i) = ourFunction(h, i); 
	}
	
	
	//x.raw_print();
	//A.print(); 

	/*v = solve(A,f);*/


	mat L, U;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	lu(L, U, A);

	vec w = solve(L, f);
	v = solve(U, w);

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	/*cout << duration_cast<chrono::milliseconds>(t1 - t2).count()<< endl;*/
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	/*v.raw_print("v=");*/

	return time_span.count(); 


	return 0; 

}
