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

void fit::non_lin_fit(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndata, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector<std::vector<double>> &covar, std::vector<std::vector<double>> &alpha, double *chisq, void(*funcs)(double, std::vector<double> &, double *, std::vector<double> &, int &), int &itmax, double &toler)
{
	// Levenberg-Marquardt method, attempting to reduce the value chi^{2} of a fit between a set of data points x[0..na-1], y[na-1]
	// with individual standard deviations sig[0..na-1] and a nonlinear function dependent on ma coefficients a[0..ma-1]
	// R. Sheehan 13 - 7 - 2018
	// Fixed the covariance matrix error and the error which meant that single-parameter fits could not be performed. 
	// Feeling pleased
	// R. Sheehan 22 - 11 - 2019

	try {
		bool c1 = ndata > 3 ? true : false;
		bool c2 = (int)(x.size()) == ndata ? true : false;
		bool c3 = (int)(y.size()) == ndata ? true : false;
		bool c4 = (int)(sig.size()) == ndata ? true : false;
		bool c5 = ma > 0 ? true : false;
		bool c6 = (int)(a.size()) == ma ? true : false;
		bool c7 = (int)(ia.size()) == ma ? true : false;
		std::pair<int, int> sze = lin_alg::array_2D_size(alpha);
		bool c9 = sze.first == ma ? true : false;
		bool c10 = sze.second == ma ? true : false;
		std::pair<int, int> szea = lin_alg::array_2D_size(covar);
		bool c9a = szea.first == ma ? true : false;
		bool c10a = szea.second == ma ? true : false;
		bool c11 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c9 && c10 && c9a && c10a ? true : false;

		if (c11) {

			int itnum = 0;
			double ochisq, alamda; 

			ochisq = *chisq; 

			while (itnum < itmax) {

				if (itnum == 0) alamda = -1.0; 

				// call mrqmin routine until convergence is achieved
				mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, chisq, funcs, &alamda); 

				// check for convergence
				if (fabs(*chisq - ochisq) < toler) { 
					
					std::cout << "Iterations converged to within tolerance after " << itnum << " steps\n"; 

					break;
				}

				ochisq = *chisq;
				
				itnum++; 
			}

			// Final call to mrqmin with alamda = 0, so that covar[1..ma][1..ma] returns the 
			// covariance matrix, and alpha the curvature matrix
			alamda = 0.0; 
			mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, chisq, funcs, &alamda);

			double nu = ndata - ma; 
			double q = probability::gammq(0.5*nu, 0.5*(*chisq)); // goodness of fit

			std::cout << "The computed fit parameters are:\n";
			int count = 0; 
			for (int i = 0; i < ma; i++) {
				if (ia[i]) {
					std::cout << "a[" << i << "] = " << a[i] << " +/- " << sqrt(covar[count][count]) << "\n";
					count++; 
				}
				else {
					std::cout << "a[" << i << "] = " << a[i] << "\n";
				}				
			}
			std::cout << "\nThe chi-sq value for the fit is " << *chisq << "\n"; 
			std::cout << "nu for the fit is " << nu << "\n"; 
			std::cout << "chi-sq / nu = "<<*chisq/nu<<"\n";
			std::cout << "goodness of fit = " << q << "\n\n";
		}
		else {
			std::string reason = "Error: fit::mrqmin()\n";
			if (!c1) reason += "No. data points is not correct ndata = " + template_funcs::toString(ndata) + "\n";
			if (!c2) reason += "x does not have correct size x.size() = " + template_funcs::toString(x.size()) + "\n";
			if (!c3) reason += "y does not have correct size y.size() = " + template_funcs::toString(y.size()) + "\n";
			if (!c4) reason += "sig does not have correct size sig.size() = " + template_funcs::toString(sig.size()) + "\n";
			if (!c5) reason += "No. fit parameters is not correct ma = " + template_funcs::toString(ma) + "\n";
			if (!c6) reason += "a does not have correct size a.size() = " + template_funcs::toString(a.size()) + "\n";
			if (!c7) reason += "ia does not have correct size ia.size() = " + template_funcs::toString(ia.size()) + "\n";
			if (!c9 || !c10) reason += "alpha does not have correct size alpha.size() = ( " + template_funcs::toString(sze.first) + " , " + template_funcs::toString(sze.second) + " )\n";
			if (!c9a || !c10a) reason += "covar does not have correct size covar.size() = ( " + template_funcs::toString(szea.first) + " , " + template_funcs::toString(szea.second) + " )\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void fit::mrqmin(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndata, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector<std::vector<double>> &covar, std::vector<std::vector<double>> &alpha, double *chisq, void(*funcs)(double, std::vector<double> &, double *, std::vector<double> &, int &), double *alamda)
{
	// Levenberg-Marquardt method, attempting to reduce the value chi^{2} of a fit between a set of data points x[0..na-1], y[na-1]
	// with individual standard deviations sig[0..na-1] and a nonlinear function dependent on ma coefficients a[0..ma-1]

	// The input array ia[0..ma-1] indicates by nonzero entries those components of a that should be fitted for, and by zero
	// entries those components that should be held fixed at their input values. 
	
	// The program returns current best - fit values for the parameters in a[0..ma-1] and chi^{2}.  
	// The arrays covar[0..ma-1][0..ma-1] and alpha[0..ma-1][0..ma-1] are used as working space during most iterations. 
	
	// Routine funcs evaluates the fitting function and its derivatives dyda[0..ma-1] with respect to the fitting parameters a at x. 
	
	// On the first call provide an initial guess for the parameters a, and set alamda<0 for initialization (which then sets 
	// alamda = .001). If a step succeeds chisq becomes smaller and alamda decreases by a factor of 10. If a step fails alamda 
	// grows by a factor of 10

	// You must call this routine repeatedly until convergence is achieved. Then, make one final 
	// call with alamda = 0, so that covar[1..ma][1..ma] returns the covariance matrix, and alpha the curvature matrix. 
	// (Parameters held fixed will return zero covariances.)

	// R. Sheehan 13 - 7 - 2018

	try {
	
		bool c1 = ndata > 3 ? true : false;
		bool c2 = (int)(x.size()) == ndata ? true : false;
		bool c3 = (int)(y.size()) == ndata ? true : false;
		bool c4 = (int)(sig.size()) == ndata ? true : false;
		bool c5 = ma > 0 ? true : false;
		bool c6 = (int)(a.size()) == ma ? true : false;
		bool c7 = (int)(ia.size()) == ma ? true : false;
		std::pair<int, int> sze = lin_alg::array_2D_size(alpha);
		bool c9 = sze.first == ma ? true : false;
		bool c10 = sze.second == ma ? true : false;
		std::pair<int, int> szea = lin_alg::array_2D_size(covar);
		bool c9a = szea.first == ma ? true : false;
		bool c10a = szea.second == ma ? true : false;
		bool c11 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c9 && c10 && c9a && c10a ? true : false;

		if (c11) {

			int j, k, l, ncols = 1;
			static int mfit;
			static double ochisq;

			std::vector<double> da;
			std::vector<double> beta; 
			std::vector<double> atry; 

			std::vector<std::vector<double>> oneda;
			std::vector<std::vector<double>> temp;

			atry = std::vector<double>(ma, 0.0);
			beta = std::vector<double>(ma, 0.0);
			da = std::vector<double>(ma, 0.0);

			for (mfit = 0, j = 0; j < ma; j++) if (ia[j]) mfit++;

			oneda = lin_alg::array_2D(mfit, ncols);
			temp = lin_alg::array_2D(mfit, mfit);

			*alamda = 0.001;

			mrqcof(x, y, sig, ndata, a, ia, ma, alpha, beta, chisq, funcs);

			ochisq = (*chisq);

			for (j = 0; j < ma; j++) atry[j] = a[j];
			
			// Alter linearised fit matrix along its diagonal elements
			for (j = 0; j < mfit; j++) {
				for (k = 0; k < mfit; k++) covar[j][k] = alpha[j][k];
				covar[j][j] = alpha[j][j] * (1.0 + (*alamda));
				for (k = 0; k < mfit; k++) temp[j][k] = covar[j][k]; 
				oneda[j][0] = beta[j];
			}

			// solve the system of equations covar . x = oneda, store solution in oneda
			// do calculation for covar using array temp
			lin_alg::gaussj(temp, mfit, oneda, ncols);
			
			// update the vector da
			for (j = 0; j < mfit; j++) {
				for (k = 0; k < mfit; k++) covar[j][k] = temp[j][k]; // retrieve computed values for covar from temp
				da[j] = oneda[j][0];
			}
			
			// if solution is converged evaluate covariance matrix
			if (*alamda == 0.0) {
				covsrt(covar, ma, ia, mfit);
				covsrt(alpha, ma, ia, mfit);

				oneda.clear(); temp.clear(); da.clear(); beta.clear(); atry.clear(); 
				
				return;
			}

			// Examine whether or not the trial is successful
			for (j = 0, l = 0; l < ma; l++)	if (ia[l]) atry[l] = a[l] + da[j++];

			// mrqcof overwrites covar, so work with temp at this point instead of covar
			// that way the correct values for covar are stored in covar
			mrqcof(x, y, sig, ndata, atry, ia, ma, temp, da, chisq, funcs);

			// If success, accept the new solution
			if (*chisq < ochisq) {
				
				*alamda *= 0.1;

				ochisq = (*chisq);

				for (j = 0; j < mfit; j++) {
					for(k=0; k < mfit; k++) alpha[j][k] = covar[j][k];
					beta[j] = da[j];
				}

				for (l = 0; l < ma; l++) a[l] = atry[l]; 
			}
			else { // solution is not accurate repeat process with updated alamda value
				*alamda *= 10.0;
				*chisq = ochisq;
			}

			oneda.clear(); temp.clear(); da.clear(); beta.clear(); atry.clear();
		}
		else {
			std::string reason = "Error: fit::mrqmin()\n";
			if (!c1) reason += "No. data points is not correct ndata = " + template_funcs::toString(ndata) + "\n";
			if (!c2) reason += "x does not have correct size x.size() = " + template_funcs::toString(x.size()) + "\n";
			if (!c3) reason += "y does not have correct size y.size() = " + template_funcs::toString(y.size()) + "\n";
			if (!c4) reason += "sig does not have correct size sig.size() = " + template_funcs::toString(sig.size()) + "\n";
			if (!c5) reason += "No. fit parameters is not correct ma = " + template_funcs::toString(ma) + "\n";
			if (!c6) reason += "a does not have correct size a.size() = " + template_funcs::toString(a.size()) + "\n";
			if (!c7) reason += "ia does not have correct size ia.size() = " + template_funcs::toString(ia.size()) + "\n";
			if (!c9 || !c10) reason += "alpha does not have correct size alpha.size() = ( " + template_funcs::toString(sze.first) + " , " + template_funcs::toString(sze.second) + " )\n";
			if (!c9a || !c10a) reason += "covar does not have correct size covar.size() = ( " + template_funcs::toString(szea.first) + " , " + template_funcs::toString(szea.second) + " )\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void fit::mrqcof(std::vector<double> &x, std::vector<double> &y, std::vector<double> &sig, int &ndata, std::vector<double> &a, std::vector<int> &ia, int &ma, std::vector<std::vector<double>> &alpha, std::vector<double> &beta, double *chisq, void(*funcs)(double, std::vector<double> &, double *, std::vector<double> &, int&))
{
	// Used by mrqmin to evaluate the linearised alpha and beta arrays and to compute chi^{2} of a fit between a set of data points 
	// x[0..na-1], y[na-1] with individual standard deviations sig[0..na-1] and a nonlinear function dependent on ma coefficients 
	// a[0..ma-1]

	// The input array ia[0..ma-1] indicates by nonzero entries those components of a that should be fitted for, and by zero
	// entries those components that should be held fixed at their input values. 

	// This functions computes the system of linear equations that must be solved for at each step
	// alpha . delta = beta

	// Routine funcs evaluates the fitting function and its derivatives dyda[0..ma-1] with respect to the fitting parameters a at x. 

	// R. Sheehan 13 - 7 - 2018

	try {

		bool c1 = ndata > 3 ? true : false;
		bool c2 = (int)(x.size()) == ndata ? true : false;
		bool c3 = (int)(y.size()) == ndata ? true : false;
		bool c4 = (int)(sig.size()) == ndata ? true : false;
		bool c5 = ma > 0 ? true : false;
		bool c6 = (int)(a.size()) == ma ? true : false;
		bool c7 = (int)(ia.size()) == ma ? true : false;
		bool c8 = (int)(beta.size()) == ma ? true : false;
		bool c11 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 ? true : false;

		if (c11) {

			int i, j, k, l, m, mfit = 0;
			double ymod, wt, sig2i, dy;

			std::vector<double> dyda(ma, 0.0);

			for (j = 0; j < ma; j++) if (ia[j]) mfit++;

			// initialise symmetric alpha, beta	
			for (j = 0; j < mfit; j++) {
				for (k = 0; k <= j; k++) alpha[j][k] = 0.0;
				beta[j] = 0.0;
			}

			*chisq = 0.0;

			// sum loop over all data
			for (i = 0; i < ndata; i++) {

				(*funcs)(x[i], a, &ymod, dyda, ma);

				sig2i = 1.0 / ( template_funcs::DSQR(sig[i]) );

				dy = y[i] - ymod;

				for (j = 0, l = 0; l < ma; l++) {
					if (ia[l]) {
						wt = dyda[l] * sig2i;
						for (k = 0, m = 0; m <= l; m++)
							if (ia[m]) alpha[j][k++] += wt * dyda[m];
						beta[j++] += dy * wt;
					}
				}

				*chisq += template_funcs::DSQR(dy)*sig2i; // find chi^{2}
			}

			// fill the symmetric side
			for (j = 1; j < mfit; j++)
				for (k = 0; k<j; k++) alpha[k][j] = alpha[j][k];

			dyda.clear();
		}
		else {
			std::string reason = "Error: fit::mrqcof()\n";
			if (!c1) reason += "No. data points is not correct ndata = " + template_funcs::toString(ndata) + "\n";
			if (!c2) reason += "x does not have correct size x.size() = " + template_funcs::toString(x.size()) + "\n";
			if (!c3) reason += "y does not have correct size y.size() = " + template_funcs::toString(y.size()) + "\n";
			if (!c4) reason += "sig does not have correct size sig.size() = " + template_funcs::toString(sig.size()) + "\n";
			if (!c5) reason += "No. fit parameters is not correct ma = " + template_funcs::toString(ma) + "\n";
			if (!c6) reason += "a does not have correct size a.size() = " + template_funcs::toString(a.size()) + "\n";
			if (!c7) reason += "ia does not have correct size ia.size() = " + template_funcs::toString(ia.size()) + "\n";
			if (!c8) reason += "beta does not have correct size beta.size() = " + template_funcs::toString(beta.size()) + "\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

void fit::goodness_of_fit()
{
	// To characterise a fit as good
	// nu = no. data points - no. fit parameters
	// chi^{2} / nu ~ 1
	// q = gammq(nu/2, chi^{2}/2) >= 0.1
	// R^{2} ~ 1
}