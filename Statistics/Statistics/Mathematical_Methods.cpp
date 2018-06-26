#ifndef ATTACH_H
#include "Attach.h"
#endif

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10

double math_meth::brent(double ax, double bx, double cx, double(*f)(double), double tol, double *xmin)
{
	// Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
	// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
	// the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
	// the minimum is returned as xmin, and the minimum function value is returned as brent, the
	// returned function value.

	try {
		int iter;
		double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm, e = 0.0;

		a = (ax < cx ? ax : cx);
		b = (ax > cx ? ax : cx);
		x = w = v = bx;
		fw = fv = fx = (*f)(x);
		for (iter = 1; iter <= ITMAX; iter++) {
			xm = 0.5*(a + b);
			tol2 = 2.0*(tol1 = tol * fabs(x) + ZEPS);
			if (fabs(x - xm) <= (tol2 - 0.5*(b - a))) {
				*xmin = x;
				return fx;
			}
			if (fabs(e) > tol1) {
				r = (x - w)*(fx - fv);
				q = (x - v)*(fx - fw);
				p = (x - v)*q - (x - w)*r;
				q = 2.0*(q - r);
				if (q > 0.0) p = -p;
				q = fabs(q);
				etemp = e;
				e = d;
				if (fabs(p) >= fabs(0.5*q*etemp) || p <= q * (a - x) || p >= q * (b - x))
					d = CGOLD * (e = (x >= xm ? a - x : b - x));
				else {
					d = p / q;
					u = x + d;
					if (u - a < tol2 || b - u < tol2)
						d = template_funcs::SIGN(tol1, xm - x);
				}
			}
			else {
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			}
			u = (fabs(d) >= tol1 ? x + d : x + template_funcs::SIGN(tol1, d));
			fu = (*f)(u);
			if (fu <= fx) {
				if (u >= x) a = x; else b = x;
				template_funcs::SHFT(v, w, x, u);
				template_funcs::SHFT(fv, fw, fx, fu);
			}
			else {
				if (u < x) a = u; else b = u;
				if (fu <= fw || w == x) {
					v = w;
					w = u;
					fv = fw;
					fw = fu;
				}
				else if (fu <= fv || v == x || v == w) {
					v = u;
					fv = fu;
				}
			}
		}
		std::string reason = "Error: math_meth::brent()\nToo many iterations\n"; 
		throw std::runtime_error(reason); 
		*xmin = x;
		return fx;
	}
	catch (std::runtime_error &e) {
		return 0.0; 
		std::cerr << e.what();
	}
}

#undef ITMAX
#undef CGOLD
#undef ZEPS

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20

void math_meth::mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double(*func)(double))
{
	// Given a function func, and given distinct initial points ax and bx, this routine searches in
	 // the downhill direction(defined by the function as evaluated at the initial points) and returns
	// new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
	// values at the three points

	try{
		bool c1 = fabs(*ax - *bx) > TINY ? true : false; 
		bool c2 = fabs(*bx - *cx) > TINY ? true : false;

		if (c1 && c2) {
			double ulim, u, r, q, fu, dum, ar4;

			*fa = (*func)(*ax);
			*fb = (*func)(*bx);
			if (*fb > *fa) {
				template_funcs::SHFT(dum, *ax, *bx, dum);
				template_funcs::SHFT(dum, *fb, *fa, dum);
			}
			*cx = (*bx) + GOLD * (*bx - *ax);
			*fc = (*func)(*cx);
			while (*fb > *fc) {
				r = (*bx - *ax)*(*fb - *fc);
				q = (*bx - *cx)*(*fb - *fa);
				u = (*bx) - ((*bx - *cx)*q - (*bx - *ax)*r) /
					(2.0*template_funcs::SIGN(std::max(fabs(q - r), TINY), q - r));
				ulim = (*bx) + GLIMIT * (*cx - *bx);
				if ((*bx - u)*(u - *cx) > 0.0) {
					fu = (*func)(u);
					if (fu < *fc) {
						*ax = (*bx);
						*bx = u;
						*fa = (*fb);
						*fb = fu;
						return;
					}
					else if (fu > *fb) {
						*cx = u;
						*fc = fu;
						return;
					}
					u = (*cx) + GOLD * (*cx - *bx);
					fu = (*func)(u);
				}
				else if ((*cx - u)*(u - ulim) > 0.0) {
					fu = (*func)(u);
					if (fu < *fc) {
						ar4 = *cx + GOLD * (*cx - *bx);
						template_funcs::SHFT(*bx, *cx, u, ar4);

						ar4 = (*func)(u);
						template_funcs::SHFT(*fb, *fc, fu, ar4);
					}
				}
				else if ((u - ulim)*(ulim - *cx) >= 0.0) {
					u = ulim;
					fu = (*func)(u);
				}
				else {
					u = (*cx) + GOLD * (*cx - *bx);
					fu = (*func)(u);
				}
				template_funcs::SHFT(*ax, *bx, *cx, u);
				template_funcs::SHFT(*fa, *fb, *fc, fu);
			}		
		}
		else {
			std::string reason = "Error: math_meth::mnbrak()\nInput bracket values not distinct\n";
			throw std::runtime_error(reason);
		}
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY

#define ITMAX 100
#define EPS 3.0e-8

double math_meth::zbrent(double(*func)(double), double x1, double x2, double tol)
{
	// Using Brent’s method, find the root of a function func known to lie between x1 and x2. The
	// root, returned as zbrent, will be refined until its accuracy is tol.

	try {
		int iter;
		double a = x1, b = x2, c = x2, d, e, min1, min2;
		double fa = (*func)(a), fb = (*func)(b), fc, p, q, r, s, tol1, xm;

		if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
			std::string reason = "Error: math_meth::zbrent()\nRoot must be bracketed in zbrent\n"; 
			throw std::runtime_error(reason); 
		}
		fc = fb;
		for (iter = 1; iter <= ITMAX; iter++) {
			if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
				c = a;
				fc = fa;
				e = d = b - a;
			}
			if (fabs(fc) < fabs(fb)) {
				a = b;
				b = c;
				c = a;
				fa = fb;
				fb = fc;
				fc = fa;
			}
			tol1 = 2.0*EPS*fabs(b) + 0.5*tol;
			xm = 0.5*(c - b);
			if (fabs(xm) <= tol1 || fb == 0.0) return b;
			if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
				s = fb / fa;
				if (a == c) {
					p = 2.0*xm*s;
					q = 1.0 - s;
				}
				else {
					q = fa / fc;
					r = fb / fc;
					p = s * (2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
					q = (q - 1.0)*(r - 1.0)*(s - 1.0);
				}
				if (p > 0.0) q = -q;
				p = fabs(p);
				min1 = 3.0*xm*q - fabs(tol1*q);
				min2 = fabs(e*q);
				if (2.0*p < (min1 < min2 ? min1 : min2)) {
					e = d;
					d = p / q;
				}
				else {
					d = xm;
					e = d;
				}
			}
			else {
				d = xm;
				e = d;
			}
			a = b;
			fa = fb;
			if (fabs(d) > tol1)
				b += d;
			else
				b += template_funcs::SIGN(tol1, xm);
			fb = (*func)(b);
		}
		std::string reason = "Error: math_meth::zbrent()\nMaximum number of iterations exceeded in zbrent\n";
		throw std::runtime_error(reason);
		return 0.0;
	}
	catch (std::runtime_error &e) {
		return 0.0;
		std::cerr << e.what();
	}
}

#undef ITMAX
#undef EPS