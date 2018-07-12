#ifndef ATTACH_H
#include "Attach.h"
#endif

void fit::lin_fit(std::vector<double> &x, std::vector<double> &y, int &ndata, std::vector<double> &sig, int &mwt, double *a, double *b, double *siga, double *sigb, double *chi2, double *q)
{
	// Given a set of data points x[1..ndata], y[1..ndata] with individual standard deviations  sig[1..ndata], 
	// fit them to a straight line y = a + bx by minimizing \chi^{2}. Returned are a, b and their respective
	// probable uncertainties siga and sigb, the chi - square chi2, and the goodness - of - fit probability q
	// (that the fit would have \chi^{2} this large or larger). If mwt=0 on input, then the standard deviations 
	// are assumed to be unavailable : q is returned as 1.0 and the normalization of chi2 is to unit standard 
	// deviation on all points.

	try {
		bool c1 = ndata > 3 ? true : false; 
		bool c2 = (int)(x.size()) > 3 ? true : false; 
		bool c3 = (int)(y.size()) > 3 ? true : false;
		bool c4 = (int)(x.size()) == ndata ? true : false;
		bool c5 = (int)(y.size()) == (int)(x.size()) ? true : false; 
		bool c6 = !mwt || mwt && (int)(sig.size()) == (int)(x.size()) ? true : false; 
		bool c7 = c1 && c2 && c3 && c4 && c5 && c6 ? true : false; 

		if (c7) {
			int i;
			double wt, t, sxoss, sx = 0.0, sy = 0.0, st2 = 0.0, ss, sigdat;

			*b = 0.0;
			if (mwt) { // compute sums with weights
				ss = 0.0;
				for (i = 0; i < ndata; i++) {
					wt = 1.0 / template_funcs::DSQR(sig[i]);
					ss += wt;
					sx += x[i] * wt;
					sy += y[i] * wt;
				}
			}
			else { // compute sums without weights
				for (i = 0; i < ndata; i++) {
					sx += x[i];
					sy += y[i];
				}
				ss = ndata;
			}

			sxoss = sx / ss;
			
			if (mwt) { // build system of equations with weights
				for (i = 0; i < ndata; i++) {
					t = (x[i] - sxoss) / sig[i];
					st2 += t * t;
					*b += t * y[i] / sig[i];
				}
			}
			else { // build system of equations without weights
				for (i = 0; i < ndata; i++) {
					t = x[i] - sxoss;
					st2 += t * t;
					*b += t * y[i];
				}
			}

			// Solve for a, b, sig_{a}, sig_{b}
			*b /= st2;
			*a = (sy - sx * (*b)) / ss;
			*siga = sqrt((1.0 + sx * sx / (ss*st2)) / ss);
			*sigb = sqrt(1.0 / st2);

			// Compute \chi^{2}
			*chi2 = 0.0;
			if (mwt == 0) { // For unweighted data evalluate typical sig using \chi^{2} and adjust std. dev. 
				for (i = 0; i < ndata; i++) {
					*chi2 += template_funcs::DSQR(y[i] - (*a) - (*b)*x[i]);
				}
				*q = 1.0;
				sigdat = sqrt((*chi2) / (ndata - 2));
				*siga *= sigdat;
				*sigb *= sigdat;
			}
			else {
				for (i = 0; i < ndata; i++) {
					*chi2 += template_funcs::DSQR((y[i] - (*a) - (*b)*x[i]) / sig[i]);
				}
				*q = probability::gammq(0.5*(ndata - 2), 0.5*(*chi2));
			}
		}
		else {
			std::string reason;
			reason = "Error: fit::lin_fit()\n";
			if (!c1) reason += "No. data points is insufficient: ndata = " + template_funcs::toString(ndata) + "\n";
			if (!c2) reason += "x.size() = " + template_funcs::toString(x.size()) + " is too small\n";
			if (!c3) reason += "y.size() = " + template_funcs::toString(y.size()) + " is too small\n";
			if (!c4) reason += "x.size() = " + template_funcs::toString(x.size()) + "and ndata = " + template_funcs::toString(ndata) + " are not equal\n";
			if (!c5) reason += "x.size() = " + template_funcs::toString(x.size()) + "and y.size() = " + template_funcs::toString(y.size()) + " are not equal\n";
			if (!c6) reason += "x.size() = " + template_funcs::toString(x.size()) + "and sig.size() = " + template_funcs::toString(y.size()) + " are not equal\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void fit::covsrt(std::vector< std::vector< double > > &covar, int &ma, std::vector<int> &ia, int &mfit)
{
	// Expand in storage the covariance matrix covar, so as to take into account parameters that are
	// being held fixed. (For the latter, return zero covariances.)

	try {

		if (lin_alg::array_2D_square(covar)) {
			int i, j, k;

			for (i = mfit; i < ma; i++)
				for (j = 0; j < i + 1; j++) covar[i][j] = covar[j][i] = 0.0;

			k = mfit - 1;
			for (j = ma - 1; j >= 0; j--) {
				if (ia[j]) {
					for (i = 0; i < ma; i++) template_funcs::SWAP(covar[i][k], covar[i][j]);
					for (i = 0; i < ma; i++) template_funcs::SWAP(covar[k][i], covar[j][i]);
					k--;
				}
			}

		}
		else {
			std::string reason = "Error: void fit::covsrt()\nInput array is not square\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fit::lfit(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndat, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector< std::vector< double > > &covar, double *chisq, void (*funcs)(double &, std::vector<double> &, int &))
{
	// Given a set of data points x[1..ndat], y[1..ndat] with individual standard deviations sig[1..ndat], use chi^{2} minimization to fit for some or all of the coefficients a[1..ma] of
	// a function that depends linearly on a, y\,=\,\sum_{i} a_{i} func_{i}(x_{i}). The input array ia[1..ma] indicates by nonzero entries those components of a that should be fitted for, 
	// and by zero entries those components that should be held fixed at their input values.The program returns values for a[1..ma], chi^{2}, and the covariance matrix covar[1..ma][1..ma]. 
	// (Parameters held fixed will return zero covariances.) The user supplies a routine funcs(x, afunc, ma) that returns the ma basis functions evaluated at x = x in the array afunc[1..ma].

	try {
	
		bool c1 = ndat > 3 ? true : false;
		bool c2 = (int)(x.size()) > 3 ? true : false;
		bool c3 = (int)(y.size()) > 3 ? true : false;
		bool c4 = (int)(x.size()) == ndat ? true : false;
		bool c5 = (int)(y.size()) == (int)(x.size()) ? true : false;
		bool c6 = (int)(sig.size()) == (int)(x.size()) ? true : false;
		bool c7 = c1 && c2 && c3 && c4 && c5 && c6 ? true : false;

		if (c7) {		 
			// test to see if there are any parameters to be fitted
			int j, mfit = 0; 
			for (j = 0; j < ma; j++)
				if (ia[j]) mfit++;

			if (mfit == 0) {
				std::string reason = "Error: void fit::lfit()\nno parameters to be fitted\n"; 
				throw std::runtime_error(reason); 
			}
			else {
				int i, k, l, m, ncols = 1;
				double ym, wt, sum, sig2i;

				std::vector<double> afunc(ma, 0.0);
				std::vector< std::vector< double > > beta;

				beta = lin_alg::array_2D(ma, ncols);

				for (j = 0; j < mfit; j++) {
					for (k = 0; k < mfit; k++) covar[j][k] = 0.0;
					beta[j][0] = 0.0;
				}

				for (i = 0; i < ndat; i++) {
					(*funcs)(x[i], afunc, ma);
					ym = y[i];
					if (mfit < ma) {
						for (j = 0; j < ma; j++)
							if (!ia[j]) ym -= a[j] * afunc[j];
					}
					sig2i = 1.0 / template_funcs::DSQR(sig[i]);
					for (j = 0, l = 0; l < ma; l++) {
						if (ia[l]) {
							wt = afunc[l] * sig2i;
							for (k = 0, m = 0; m <= l; m++)
								if (ia[m]) covar[j][k++] += wt * afunc[m];
							beta[j++][0] += ym * wt;
						}
					}
				}

				for (j = 1; j < mfit; j++)
					for (k = 0; k<j; k++)
						covar[k][j] = covar[j][k];

				lin_alg::gaussj(covar, mfit, beta, ncols);

				for (j = 0, l = 0; l < ma; l++)
					if (ia[l]) a[l] = beta[j++][0];

				*chisq = 0.0;
				for (i = 0; i < ndat; i++) {
					(*funcs)(x[i], afunc, ma);
					for (sum = 0.0, j = 0; j < ma; j++) sum += a[j] * afunc[j];
					*chisq += template_funcs::DSQR((y[i] - sum) / sig[i]);
				}

				covsrt(covar, ma, ia, mfit);

				afunc.clear(); beta.clear();				
			}
		}
		else {
			std::string reason;
			reason = "Error: fit::lin_fit()\n";
			if (!c1) reason += "No. data points is insufficient: ndata = " + template_funcs::toString(ndat) + "\n";
			if (!c2) reason += "x.size() = " + template_funcs::toString(x.size()) + " is too small\n";
			if (!c3) reason += "y.size() = " + template_funcs::toString(y.size()) + " is too small\n";
			if (!c4) reason += "x.size() = " + template_funcs::toString(x.size()) + "and ndata = " + template_funcs::toString(ndat) + " are not equal\n";
			if (!c5) reason += "x.size() = " + template_funcs::toString(x.size()) + "and y.size() = " + template_funcs::toString(y.size()) + " are not equal\n";
			if (!c6) reason += "x.size() = " + template_funcs::toString(x.size()) + "and sig.size() = " + template_funcs::toString(y.size()) + " are not equal\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
	}
}

void fit::mrqcof(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndata, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector<std::vector<double>> &alpha, std::vector<double> &beta, double *chisq, void(*funcs)(double, std::vector<double> &, double *, std::vector<double> &, int&))
{
	// Used by mrqmin to evaluate the linearised alpha and beta arrays and to compute chi^{2}

	try {
	
		bool c1 = ndata > 3 ? true : false; 
		bool c2 = (int)(x.size()) == ndata ? true : false; 
		bool c3 = (int)(y.size()) == ndata ? true : false;
		bool c4 = (int)(sig.size()) == ndata ? true : false;
		bool c5 = ma > 0 ? true : false;
		bool c6 = (int)(a.size()) == ma ? true : false;
		bool c7 = (int)(ia.size()) == ma ? true : false;
		bool c8 = (int)(beta.size()) == ma ? true : false;
		std::pair<int, int> sze = lin_alg::array_2D_size(alpha); 
		bool c9 = sze.first == ma ? true : false;
		bool c10 = sze.second == ma ? true : false;
		bool c11 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 && c9 && c10 ? true : false; 

		if (c11) {

			int i, j, k, l, m, mfit = 0;
			double ymod, wt, sig2i, dy;

			std::vector<double> dyda(ma, 0.0);

			for (j = 0; j < ma; j++)
				if (ia[j]) mfit++;

			// initialise symmetric alpha, beta	
			for (j = 0; j < mfit; j++) {
				for (k = 0; k <= j; k++) alpha[j][k] = 0.0;
				beta[j] = 0.0;
			}

			*chisq = 0.0;

			// sum loop over all data
			for (i = 0; i < ndata; i++) {
				(*funcs)(x[i], a, &ymod, dyda, ma);
				sig2i = 1.0 / (sig[i] * sig[i]);
				dy = y[i] - ymod;
				for (j = 0, l = 0; l < ma; l++) {
					if (ia[l]) {
						wt = dyda[l] * sig2i;
						for (j++, k = 0, m = 0; m <= l; m++)
							if (ia[m]) alpha[j][++k] += wt * dyda[m];
						beta[j] += dy * wt;
					}
				}
				*chisq += dy * dy*sig2i; // find chi^{2}
			}

			// fill the symmetric side
			for (j = 1; j < mfit; j++)
				for (k = 0; k<j; k++) alpha[k][j] = alpha[j][k];

			dyda.clear();		
		}
		else {
			std::string reason = "Error: fit::mrqcof()\n";
			
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	} 
}

void fit::mrqmin(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndata, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector<std::vector<double>> &covar, std::vector<std::vector<double>> &alpha, double *chisq, void(*funcs)(double, std::vector<double> &, double *, std::vector<double> &, int &), double *alamda)
{
	
}