#ifndef H_MAPP
#define H_MAPP

double Factorial(int n);
double Combination(int n, int k);
void ShowComb(int n, int k, std::string up, std::string down);
std::string NewtoneBinomeS(int n, std::string s);
double BernsteinPoly(int n, int k, double t);
std::vector<double> VectorProduct(const std::vector<double> x, const std::vector<double> y);
double VectorProductZ(const std::vector<double> x, const std::vector<double> y);
double ScalarProduct(const std::vector<double> x, const std::vector<double> y);
double ScalarProduct2(const std::vector<double> x,const std::vector<double> y);

#endif
