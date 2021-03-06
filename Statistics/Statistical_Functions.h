#ifndef STATISTICAL_FUNCTIONS_H
#define STATISTICAL_FUNCTIONS_H

// Declaration of a namespace that contains functions used to perform statistical tests on data sets
// R. Sheehan 5 - 9 - 2017

namespace statistic{
	
	// compute the average and variance of a &data set
	void avevar(std::vector<double> &data, int &n, double *ave, double *var); 

	// compute the first four moments of a &data set
	void moment(std::vector<double> &data, int &n, double *ave, double *adev, double *sdev, double *var, double *skew, double *curt); 

	// perform t-test to check for significantly different means on two &data sets with equal variance
	void ttest(std::vector<double> &data1, int &n1, std::vector<double> &data2, int &n2, double *t, double *prob); 

	// perform t-test to check for significantly different means on two &data sets with unequal variance
	void tutest(std::vector<double> &data1, int &n1, std::vector<double> &data2, int &n2, double *t, double *prob);

	// perform t-test to check for significantly different means on two &data sets which are paired
	void tptest(std::vector<double> &data1, std::vector<double> &data2, int &n, double *t, double *prob); 

	// perform f-test to check for significantly different variances on two &data sets
	void ftest(std::vector<double> &data1, int &n1, std::vector<double> &data2, int &n2, double *f, double *prob);

	// perform a chi-square test to compare the data in bins with the data in ebins
	void chsone(std::vector<double> &bins, std::vector<double> &ebins, int &nbins, int &knstrn, double *df, double *chsq, double *prob); 

	// perform a chi-square test to compare the data in bins1 with the data in bins2
	void chstwo(std::vector<double> &bins1, std::vector<double> &bins2, int &nbins, int &knstrn, double *df, double *chsq, double *prob); 

	// perform a chi-square test to compare the data in bins1 with the data in bins2, data have different distributions
	void chsthree(std::vector<double> &bins1, std::vector<double> &bins2, int &nbins1, int &nbins2, int &knstrn, double *df, double *chsq, double *prob);

	// perform Kolmorgorv-Smirnov test to determine if a given data set is distributed according to a known distribution
	void ksone(std::vector<double> &data, int &n, double (*func)(double), double *d, double *prob); 

	// perform Kolmorgorov-Smirnov test to determine if two data sets are drawn from the same distribution
	void kstwo(std::vector<double> &data1, int &n1, std::vector<double> &data2, int &n2, double *d, double *prob); 

	// perform contingency analysis using the Cramer's V and contingency coefficient measures of association
	void cntab1(std::vector<std::vector<int>> &nn, int &ni, int &nj, double *chisq, double *df, double *prob, double *cramrv, double *ccc); 
	
	// perform contingency analysis using the entropy based measures of association
	void cntab2(std::vector<std::vector<int>> &nn, int &ni, int &nj, double *h, double *hx, double *hy, double *hygx, double *hxgy, double *uygx, double *uxgy, double *uxy); 

	// perform Pearon's r calculation of linear correlation between two data sets
	void pearsn(std::vector<double> &x, std::vector<double> &y, int &n, double *r, double *prob, double *z);

	// perform Spearman rank correlation calculation for two data sets
	void spear(std::vector<double> &data1, std::vector<double> &data2, int &n, double *d, double *zd, double *probd, double *rs, double *probrs);

	void crank(int &n, std::vector<double> &w, double *s); // routine needed by spear

	// perform Kendall tau correlation calculation for two data sets
	void kendl1(std::vector<double> &data1, std::vector<double> & data2, int &n, double *tau, double *z, double *prob); 

	// perform Kendall tau correlation calculation for a contingency table
	void kendl2(std::vector<std::vector<double>> &tab, int &i, int &j, double *tau, double *z, double *prob); 
}

#endif