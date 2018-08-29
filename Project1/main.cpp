#include <iostream>
#include <cmath>
#include <fstream>
#include <istream>
#include <string>
#include <vector>
#include <ctime>
#include <chrono>

#include "Gaussian.h"


using namespace std;

int const n = 1000;


int main() {

	//int n; 
	//cout << "Type the number of grid points you want, n" << endl; 
	//cin >> n; 

	problem1b(n); 
	/*problem1c(n); */

	//const int COLUMNS = 3; 

	//vector <vector <double>> data; 
	//string filename = "matrixdata2.txt";
	//ifstream ifile(filename.c_str());

	//if (ifile.is_open()) {
	//	int num; 

	//	vector <double> numbers_in_line; 

	//	while (ifile >> num) {
	//		numbers_in_line.push_back(num); 
	//		if (numbers_in_line.size() == COLUMNS) {
	//			data.push_back(numbers_in_line); 
	//			numbers_in_line.clear(); 
	//			}
	//	}
	//}
	//else {
	//	cerr << "Error opening the file \n"; 
	//	exit(1);
	//}
	//
	//vector <double> column; 
	//int col = 1; 

	//for (int i = 0; i < data.size(); i++) {
	//	column.push_back(data[i][col - 1]); 
	//	cout << column[i] << endl; 
	//}


	return 0; 
}

