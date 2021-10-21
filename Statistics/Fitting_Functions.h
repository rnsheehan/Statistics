#ifndef FITTING_FUNCTIONS_H
#define FITTING_FUNCTIONS_H

// Namespace containing fitting functions for linear, polynomial and nonlinear fits
// R. Sheehan 22 - 6 - 2018

namespace fit {
	// Perform linear for data with errors in y measurements only
	void lin_fit(std::vector<double> &x, std::vector<double> &y, int &ndata, std::vector<double> &sig, int &mwt, double *a, double *b, double *siga, double *sigb, double *chi2, double *q); 

	// Function for expanding the storage of an array
	void covsrt(std::vector< std::vector< double > > &covar, int &ma, std::vector<int> &ia, int &mfit); 

	// Perform least squares polynomial fit to data with errors in y values
	void lfit(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndat, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector< std::vector< double > > &covar, double *chisq, void(*funcs)(double &, std::vector<double> &, int&)); 

	// Perform fit to non-linear function
	void non_lin_fit(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndata, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector<std::vector<double>> &covar, std::vector<std::vector<double>> &alpha, double *chisq, void(*funcs)(double, std::vector<double> &, double *, std::vector<double> &, int &), int &itmax, double &toler, bool loud = false);

	// Levenberg-Marquart Algorithm
	void mrqmin(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndata, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector<std::vector<double>> &covar, std::vector<std::vector<double>> &alpha, double *chisq, void(*funcs)(double, std::vector<double> &, double *, std::vector<double> &, int &), double *alamda);

	// Function for computing systems of equations used in LMA
	void mrqcof(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndata, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector<std::vector<double>> &alpha, std::vector<double> &beta, double *chisq, void(*funcs)(double, std::vector<double> &, double *, std::vector<double> &, int &));

	// Function for evaluating goodness-of-fit statistics
	void goodness_of_fit(); 

	// Function for computing the residual of a fit
	void residuals(std::vector<double>& x, std::vector<double>& y, std::vector<double>& sig, int& ndata, std::vector<double>& a, int& ma, void(*funcs)(double, std::vector<double>&, double*, std::vector<double>&, int&), std::vector<std::vector<double>> &data);

	// generate a randomly sampled data set from an existing data set
	void random_sample(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndata, std::vector<double> &xs, std::vector<double> &ys, std::vector<double> &sigs, int &ns);
}

#endif
