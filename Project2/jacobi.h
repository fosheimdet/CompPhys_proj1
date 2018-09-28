

arma::mat constructA(double rho_min, double rho_max, int n, bool potential, bool interacting, double omega);
double maxOffDiag(arma::mat & A, int &k, int &l, int n);
void rotate(arma::mat &A, arma::mat &R, int k, int l, int n);
void jacobiMethod(arma::mat &A, arma::mat &R, int n, arma::vec &eigVecCol);
arma::vec sortEigenvalues(arma::mat A, int n);
void writeToFile(arma::mat &R, double rho_min, double rho_max, int n, int colIndex, std::string filename);
std::tuple<double,double> jacobiVsArmadillo(arma::mat &A, arma::mat& R, arma::cx_vec eigVal, arma::cx_mat eigVec, arma::vec &eigVecCol, int n);
