#ifndef ATTACH_H
#include "Attach.h"
#endif

std::vector< std::vector< double > > lin_alg::array_2D(int &nrows, int &ncols)
{
	// create an array of zeroes of size nrows*ncols
	// R. Sheehan 22 - 6 - 2018

	try {
		if (nrows > 0 && ncols > 0) {
			std::vector< std::vector< double > > name; 
			name.resize(nrows);
			for (int i = 0; i < nrows; i++) name[i].resize(ncols,0.0); 
			return name; 
		}
		else {
			std::vector< std::vector< double > > temp;
			return temp;
			std::string reason = "Error: std::vector< std::vector< double > > lin_alg::array_2D\n"; 
			if (nrows <= 1) reason += "nrows = " + template_funcs::toString(nrows) + " too small\n"; 
			if (ncols <= 1) reason += "ncols = " + template_funcs::toString(ncols) + " too small\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

std::pair<int, int> lin_alg::array_2D_size(std::vector< std::vector< double > > &name)
{
	// return a pair structure of the form (nrows, ncols)
	// telling you the size of the array
	// R. Sheehan 22 - 6 - 2018

	try {
		if ((int)(name.size()) > 0 && (int)(name[0].size()) > 0) {
			return std::pair<int, int>((int)(name.size()), (int)(name[0].size()));
		}
		else {
			return std::pair<int, int>(-1, -1);
			std::string reason = "Error: std::pair<int, int> lin_alg::array_2D_size()\n"; 
			reason += "Input array dimensions not assigned correctly\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		return std::pair<int, int>(-1, -1);
		std::cerr << e.what();
	}
}

bool lin_alg::array_2D_square(std::vector< std::vector< double > > &name)
{
	// test the 2D array to see if it is square
	// R. Sheehan 22 - 6 - 2018

	std::pair<int, int> sizes = array_2D_size(name);
		
	return sizes.first == sizes.second ? true : false; 
}

void lin_alg::gaussj(std::vector< std::vector< double > > &a, int &n, std::vector< std::vector< double > > &b, int &m)
{
	// Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
	// is the input matrix.b[1..n][1..m] is input containing the m right - hand side vectors. On
	// output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
	// vectors

	try {
		std::pair<int, int> sizea, sizeb; 

		sizea = array_2D_size(a); sizeb = array_2D_size(b);

		bool c1 = sizea.first == sizea.second ? true : false; // test for squareness of a
		bool c2 = sizea.second == sizeb.first ? true : false; // columns of A equal to rows of B? 

		if (c1 && c2) {
			int i, icol, irow, j, k, l, ll;
			double big, dum, pivinv;

			std::vector<int> indxc(n, 0);
			std::vector<int> indxr(n, 0);
			std::vector<int> ipiv(n, 0);

			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j <	n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (fabs(a[j][k]) >= big) {
									big = fabs(a[j][k]);
									irow = j;
									icol = k;
								}
							}
							else if (ipiv[k] > 1) { 
								std::string reason = "Error: void lin_alg::gaussj()\nSingular Matrix-1\n"; 
								throw std::runtime_error(reason);
							}
						}
				++(ipiv[icol]);
				
				if (irow != icol) {
					for (l = 0; l < n; l++) template_funcs::SWAP(a[irow][l], a[icol][l]);
					for (l = 0; l < m; l++) template_funcs::SWAP(b[irow][l], b[icol][l]);
				}
				
				indxr[i] = irow;
				
				indxc[i] = icol;
				
				if (a[icol][icol] == 0.0) {
					std::string reason = "Error: void lin_alg::gaussj()\nSingular Matrix-2\n";
					throw std::runtime_error(reason);
				}
				
				pivinv = 1.0 / a[icol][icol];
				
				a[icol][icol] = 1.0;
				
				for (l = 0; l < n; l++) a[icol][l] *= pivinv;
				
				for (l = 0; l < m; l++) b[icol][l] *= pivinv;
				
				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a[ll][icol];
						a[ll][icol] = 0.0;
						for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
						for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
					}
			}

			// End of Elimination procedure, unscramble row permuations
			for (l = n-1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						template_funcs::SWAP(a[k][indxr[l]], a[k][indxc[l]]);
			}
			
			ipiv.clear(); indxr.clear(); indxc.clear(); 
		}
		else {
			std::string reason = "Error: void lin_alg::gaussj()\n";
			reason += "Input array dimensions are not correct\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}