#ifndef TEST_ROUTINES_H
#define TEST_ROUTINES_H

// Declaration of namespace containing functions to test the statistical routines
// Some sample data has been provided, this is generated using Python script Data_Generator.py
// R. Sheehan 11 - 9 - 2017

namespace testing{
	
	void compute_moments(); 

	void perform_t_test_diff_mean_same_var(); 

	void perform_t_test_diff_mean_diff_var(); 

	void perform_f_test_diff_mean_diff_var(); 

	void perform_chsone_test(); 

	void perform_chstwo_test(); 

	void perform_ksone_test(); 

	void perform_kstwo_test(); 

	void perform_con_tab_test(); 	

	double ksfunc(double x); 

	void read_tab1_data(std::string &header, std::vector<std::string> &months, std::vector<std::string> &accident_types, std::vector<std::vector<int>> &numbers); 

	void read_tab2_data(std::string &header, std::vector<std::string> &months, std::vector<std::string> &cities, std::vector<std::string> &states, std::vector<double> &state_averages, std::vector<double> &state_latitudes, std::vector<std::vector<double>> &numbers); 

	void perform_pearson_test(); 

	void perform_spearman_test(); 

	void rank_order_test(); 

	void perform_kendall_test_1(); 

	void perform_kendall_test_2(); 

	void lin_fit_test(); 

	void gaussj_test(); 

	void funcs(double &x, std::vector<double> &afunc, int &ma); 

	void fpoly(double &x, std::vector<double> &afunc, int &ma);

	void fleg(double &x, std::vector<double> &afunc, int &ma);

	void lst_sqr_test(); 

	void fgauss(double x, std::vector<double> &a, double *y, std::vector<double> &dyda, int &na); 

	void non_lin_fit_test(); 

	void rng_test(); 

	void monte_carlo_test(); 
}

#endif