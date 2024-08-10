#include <iostream>
#include <valarray>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <random>

using namespace std;

double generateRandomNumber() { // Function that enables the generation of a random number between -1 and 1. This allows for the Monte Carlo method of generating multitude portfolios and checking them against
				// the criteria of low risk portfolios.
        random_device rd;
        mt19937 gen(rd());

        uniform_real_distribution<> dis(-1.0,1.0);

        return dis(gen);
}

double portfolio_variance(double w1, double w2, double v1, double v2, double v3, double c12, double c13, double c23) { // Function that calculates the variance of a constructed portfolio from three stocks.

        double var = w1*w1*v1 + w2*w2*v2 + (1 - w1 - w2)*(1 - w1 - w2) *v3 + 2*w1*w2*c12 + 2*(1 - w1 - w2)*w1*c13 + 2*w2*(1 - w1 - w2)*c23;
        return var;

}

double Lagrange1(double w1, double w2, double v1, double v3, double c12, double c13, double c23) { // Function that computes the first Lagrangian derivative, which enables the program to reach a low risk portfolio.

        return 2*w1*v1 - 2*(1 - w1 - w2)*v3 + 2*w2*c12 + (2 - 4*w1 - 2*w2)*c13 - 2*w1*c23;

}

double Lagrange2(double w1, double w2, double v2, double v3, double c12, double c13, double c23) { // Function that computes the second Lagrangian derivative, which enables the program to reach a low risk portfolio.
        return 2*w2*v2 - 2*(1 - w1 - w2)*v3 + 2*w1*c12 - 2*w1*c13 + (2 - 2*w1 - 4*w2)*c23;

}
double expected_return(valarray<double> returns) { // Function that returns the expected daily return of a given stock price.
	
	int N = returns.size();

	int i;
	double sum = 0;
	double exp;

	for (i = 0; i < N; i++) {
		sum += returns[i];

	}
	exp = sum/(double) N;
	return exp;

}

double variance_return(valarray<double> returns) { // Function that returns the variance of the daily returns of a given stock price
	int n = returns.size();
	double var = 0;
	

	int j;

	double mu = expected_return(returns);

	for (j = 0; j < n; j++) {
		var += (returns[j] - mu) * (returns[j] - mu);
	}

	double variance = var/ (double) (n-1);

	return variance;

}

double cov(valarray<double> x, valarray<double> y) { // Function that returns the covariance of the daily returns.

	if (x.size() - y.size() != 0) {
		cerr << "Error: Data sets do not have a matching number of data" << endl;
	
	} 
	int n = x.size();
	double sum = 0;

	double mux = expected_return(x);
	double muy = expected_return(y);
	int k;
	for (k = 0; k < n; k++) {
		
		sum += (x[k] - mux)*(y[k] - muy);
	}
	double covariance = sum/(double)(n-1);
	return covariance;

}


valarray<double> read_data(const string& filename, int N) { // Function that reads the given data (preferably the adjusted close performance of three stocks.
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return valarray<double>();
    }
    valarray<double> data(N * 3);
    string line;
    int row = 0;
    while (getline(file, line) && row < N) {
        stringstream ss(line);
        string cell;
        int col = 0;
        while (getline(ss, cell, ' ') && col < 3) {
            try {
                data[row * 3 + col] = stod(cell);
            } catch (const invalid_argument& e) {
                cerr << "Error: Invalid data in file " << filename << " at row " << row << endl;
                return valarray<double>();
            }
            ++col;
        }
        ++row;
    }
    return data;
}



int main(void) {

	int N = 253; // Number of days with a given adjusted close [This can be changed according to the time period of the analysed stocks]
	valarray<double> data = read_data("original.dat",N); // Data being read and processed for the program.
	valarray<double> x(N), y(N), z(N); // Value arrays for adjusted close performances of the three stocks.

	int j;

	for (j = 0; j < N; j++) { // For-loop that allocates the different adjusted closed prices of the stock for every day of the data.
		x[j] = data[3*j];
		y[j] = data[3*j + 1];
		z[j] = data[3*j + 2];

	}

	valarray<double> return_x(N-1), return_y(N-1), return_z(N-1); // Value arrays for the daily returns of each stock.

	int k;

	for (k = 0; k < N-1; k++) { // For-loop that allocates the different daily returns for the three stocks.

		return_x[k] = log(x[k+1]) - log(x[k]);
		return_y[k] = log(y[k+1]) - log(y[k]);
		return_z[k] = log(z[k+1]) - log(z[k]);

	}

	double e1, e2, e3, c12, c13, c23, v1, v2, v3; // All the statistical attributes for the three stocks: expected returns, covariance of returns, variance of returns.

	e1 = expected_return(return_x);
	e2 = expected_return(return_y);
	e3 = expected_return(return_z);

	c12 = cov(return_x, return_y);
	c13 = cov(return_x, return_z);
	c23 = cov(return_y, return_z);

	v1 = variance_return(return_x);
	v2 = variance_return(return_y);
	v3 = variance_return(return_z);
	
	int Ndata = 100000; // Number of iterations for the Monte Carlo method
        

        double TOL = 0.00001; // Tolerance for the while loop which will stop whenever both derivatives are very close to 0.
        valarray<double> w1(Ndata), w2(Ndata);

        int f;

        for (f = 0; f < Ndata; f++) { // For-loop that generates a random number between -1 and 1, which will constitute the weight (or the proportion of wealth) that shall be invested in a stock for creating 
				      // a low risk portfolio.
                w1[f] = generateRandomNumber();
                w2[f] = generateRandomNumber();

        }



        // Time to use the code to find the minimum variance portfolio

        int i = -1;


        do { // Do-while loop that will check every single iteration of possible portfolio until the Lagrangian derivatives reach a number very close to zero, which will mean that a low risk portfolio has been
	     // found for the given stocks.
                i++;
                cerr << "w1 : " << w1[i] << "  " << "w2 : " << w2[i] << "  " << "Variance : " << portfolio_variance(w1[i],w2[i],v1,v2,v3,c12,c13,c23) << " " << "DVar(Rp)/Dw1 = " << Lagrange1(w1[i],w2[i],v1,v3,c12,c13,c23) << endl;
		cerr << "DVar(Rp)/Dw2 = " << Lagrange2(w1[i], w2[i], v2, v3, c12, c13, c23) << endl;



        } while (((abs(Lagrange1(w1[i],w2[i],v1,v3,c12,c13,c23)) > TOL || abs(Lagrange2(w1[i],w2[i],v2,v3,c12,c13,c23)) > TOL)) || (e1*w1[i] + e2*w2[i] + e3*(1 - w1[i] - w2[i]) < 0) && i <= Ndata);
        cerr << "Low risk portfolio consists of : w1 = " << w1[i] << " and w2 = " << w2[i] << "   with expectation of portfolio : " << e1*w1[i] + e2*w2[i] + e3*(1-w1[i] - w2[i]) << endl;
	cerr << "With portfolio variance: " << portfolio_variance(w1[i],w2[i],v1,v2,v3,c12,c13,c23) << endl;

	cerr << "Ended after " << i+1 << " iterations." << endl;
	cerr << "Expectation of Stock 1: " << e1 << " Variance of Stock 1: " << v1 << endl;
	cerr << "Expectation of Stock 2: " << e2 << " Variance of Stock 2: " << v2 << endl;
	cerr << "Expectation of Stock 3: " << e3 << " Variance of Stock 3: " << v3 << endl;
	


return 0;



}
