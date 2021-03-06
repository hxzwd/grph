#include "main.hpp"
#include "mapp.hpp"

double Factorial(int n)
{
	double result = 1.0;
	if(n == 0 || n == 1) return result;
	for(int i = 1; i <= n; i++)
		result = result * i;
	return result;
}

double Combination(int n, int k)
{
	return Factorial(n)/(Factorial(k)*Factorial(n - k));
}

void ShowComb(int n, int k, std::string up = "", std::string down = "")
{
	if(up == "r" && down == "r") 
	{ 
		up = down = ""; 
		std::cout << "Calculated: " << Combination(n, k) << "\n"; 
	}
	std::cout << "(" << n << ") " << up << "\n";
	std::cout << "(" << k << ") " << down << "\n";
}

std::string NewtoneBinomeS(int n, std::string s = "1 + x")
{
	std::string first = "1";
	std::string second = "x";
	
	std::string coeff = "";
	std::string variable = "";
	std::string op = "+";
	std::string result = "";
	
	if(n == 0) return "1";
	if(n == 1) return s;
	
	for(int i = 0; i <= n; i++)
	{
		coeff = std::to_string((int)Combination(n, i));
		variable = "x^" + std::to_string(i);
		if(i == 1) variable = "x";
		if(i == 0) { result = result + coeff + std::string(" ") + op + std::string(" "); continue; }
		if(i == n) op = "";
		result = result + coeff + std::string("*") + variable + std::string(" ") + op + std::string(" ");
	}
	return result;
}

double BernsteinPoly(int n, int k, double t)
{	
		return Combination(n, k) * std::pow(1.0f - t, n - k) * std::pow(t, k);
}

std::vector<double> VectorProduct(std::vector<double> x, std::vector<double> y)
{
	std::vector<double> result;
	std::vector<double> X = x;
	std::vector<double> Y = y;
	if(x.size() < 3 || y.size() < 3)
	{
		X.push_back(0);
		Y.push_back(0);
	}
	
	result.push_back(X[1] * Y[2] - Y[2] * Y[1]);
	result.push_back(-(X[0] * Y[2] - Y[0] * X[2]));
	result.push_back(X[0] * Y[1] - Y[0] * X[1]);
	return result;
}

double ScalarProduct(const std::vector<double> x,const std::vector<double> y)
{
	double result = x[0] * y[0] + x[1] * y[1];
	if(x.size() == y.size() == 3)
		result = result + x[2] * y[2];
	return result;
}

double ScalarProduct2(const std::vector<double> x,const std::vector<double> y)
{
	return x[0] * y[0] + x[1] * y[1];
}

double VectorProductZ(const std::vector<double> x,const std::vector<double> y)
{
	std::vector<double> result;
	std::vector<double> X = x;
	std::vector<double> Y = y;
	if(x.size() == 2 || y.size() == 2)
	{
		X.push_back(0);
		Y.push_back(0);
	}
	result.push_back(X[1] * Y[2] - Y[2] * Y[1]);
	result.push_back(-(X[0] * Y[2] - Y[0] * X[2]));
	result.push_back(X[0] * Y[1] - Y[0] * X[1]);
	return result[2];
}
