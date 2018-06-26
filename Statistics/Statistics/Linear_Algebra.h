#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

// Implementation of functions required for linear algebra calculations inside the statistical functions
// R. Sheehan 22 - 6 - 2018

namespace lin_alg {
	
	bool array_2D_square(std::vector< std::vector< double > > &name); // test to see if an array is square

	std::vector< std::vector< double > > array_2D(int &nrows, int &ncols); // create a 2D array of given size
	
	std::pair<int, int> array_2D_size(std::vector< std::vector< double > > &name); // return the size of an array as a pair

	// Implementation of Gaussian elimination with partial pivoting 
	void gaussj(std::vector< std::vector< double > > &a, int &n, std::vector< std::vector< double > > &b, int &m);
}

#endif