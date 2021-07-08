#ifndef MATHEMATICAL_METHODS_H
#define MATHEMATICAL_METHODS_H

// Namespace containing mathematical methods for performing calculations inside the statistical routines
// R. Sheehan 22 - 6 - 2018

namespace math_meth {

	double brent(double ax, double bx, double cx, double(*f)(double), double tol, double *xmin); // Brent's minima search method

	double zbrent(double(*func)(double), double x1, double x2, double tol); // Brent's root finding method
	
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double(*func)(double)); // Minimum bracketing 
}

#endif
