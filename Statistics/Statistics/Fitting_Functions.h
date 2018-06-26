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
}

#endif
