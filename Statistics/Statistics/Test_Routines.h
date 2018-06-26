#ifndef TEST_ROUTINES_H
#define TEST_ROUTINES_H

// Declaration of namespace containing functions to test the statistical routines
// Some sample data has been provided, this is generated using Python script Data_Generator.py
// R. Sheehan 11 - 9 - 2017

namespace testing{
	
	void compute_moments(); 

	void t_test_diff_mean_same_var(); 

	void t_test_diff_mean_diff_var(); 

	void f_test_diff_mean_diff_var(); 

	void chsone_test(); 

	void chstwo_test(); 

	void ksone_test(); 

	void kstwo_test(); 

	void con_tab_test(); 	

	double ksfunc(double x); 

	void read_tab1_data(std::string &header, std::vector<std::string> &months, std::vector<std::string> &accident_types, std::vector<std::vector<int>> &numbers); 

	void read_tab2_data(std::string &header, std::vector<std::string> &months, std::vector<std::string> &cities, std::vector<std::string> &states, std::vector<double> &state_averages, std::vector<double> &state_latitudes, std::vector<std::vector<double>> &numbers); 

	void pearson_test(); 

	void spearman_test(); 

	void rank_order_test(); 

	void kendall_test_1(); 

	void kendall_test_2(); 

	void lin_fit_test(); 

	void gaussj_test(); 

	void funcs(double &x, std::vector<double> &afunc, int &ma); 

	void fpoly(double &x, std::vector<double> &afunc, int &ma);

	void fleg(double &x, std::vector<double> &afunc, int &ma);

	void lst_sqr_test(); 
}

#endif