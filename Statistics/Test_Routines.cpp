#ifndef ATTACH_H
#include "Attach.h"
#endif

// info on the pragma directive
// https://docs.microsoft.com/en-us/cpp/preprocessor/pragma-directives-and-the-pragma-keyword?view=msvc-160
// https://gcc.gnu.org/onlinedocs/cpp/Pragmas.html
#pragma warning(disable : 4834) // turn off warning that tells you to use the return value of std::remove_if
#pragma warning(disable : 26451) // turn off warning that tells you about possible arithmetic overflow


void testing::compute_moments()
{
	// Compute the moments of some of the sample data sets
	// R. Sheehan 11 - 9 - 2017

	std::string filename; 
	std::vector<double> data;
	int n_data; 
	double mean, var;

	filename = "Data_1.txt";  
	
	useful_funcs::read_into_vector(filename, data, n_data, true); 

	statistic::avevar(data, n_data, &mean, &var);

	std::cout<<"\nFile:\t" + filename + "\n"; 
	std::cout<<"Size:\t" + template_funcs::toString(n_data) + "\n"; 
	std::cout<<"Average:\t" + template_funcs::toString(mean, 4) + "\n"; 
	std::cout<<"Variance:\t" + template_funcs::toString(var, 4) + "\n\n"; 

	data.clear(); 

	filename = "Data_6.txt"; 
	
	useful_funcs::read_into_vector(filename, data, n_data, true); 

	statistic::avevar(data, n_data, &mean, &var);

	std::cout<<"\nFile:\t" + filename + "\n"; 
	std::cout<<"Size:\t" + template_funcs::toString(n_data) + "\n"; 
	std::cout<<"Average:\t" + template_funcs::toString(mean, 4) + "\n"; 
	std::cout<<"Variance:\t" + template_funcs::toString(var, 4) + "\n\n"; 

	data.clear(); 
}

void testing::perform_t_test_diff_mean_same_var()
{
	// perform a Student's t-test on two distributions with different means but the same variance
	// determine if the data have significantly different mean
	// R. Sheehan 11 - 9 - 2017

	std::string file1, file2; 
	std::vector<double> data1, data2;
	int n_data_1, n_data_2;
	double mean, var, t_stat, p_val; 

	bool loud = true; 

	//file1 = "Data_1.txt"; // mu = 0, sig = 1
	//file2 = "Data_6.txt";  // mu = 1, sig = 1

	//file1 = "Data_3.txt"; // mu = 1.5, sig = 0.5
	file1 = "Data_4.txt";  // mu = 0, sig = 0.5
	file2 = "Data_5.txt";  // mu = 0.01, sig = 0.5

	useful_funcs::read_into_vector(file1, data1, n_data_1, loud);

	statistic::avevar(data1, n_data_1, &mean, &var);

	std::cout<<"\nFile:\t" + file1 + "\n"; 
	std::cout<<"Size:\t" + template_funcs::toString(n_data_1) + "\n"; 
	std::cout<<"Average:\t" + template_funcs::toString(mean, 4) + "\n"; 
	std::cout<<"Variance:\t" + template_funcs::toString(var, 4) + "\n\n";

	useful_funcs::read_into_vector(file2, data2, n_data_2, loud);

	statistic::avevar(data2, n_data_2, &mean, &var);

	std::cout<<"\nFile:\t" + file2 + "\n"; 
	std::cout<<"Size:\t" + template_funcs::toString(n_data_2) + "\n"; 
	std::cout<<"Average:\t" + template_funcs::toString(mean, 4) + "\n"; 
	std::cout<<"Variance:\t" + template_funcs::toString(var, 4) + "\n\n";

	statistic::ttest(data1, n_data_1, data2, n_data_2, &t_stat, &p_val);

	std::cout<<"\nResult of Student's t-test, assuming equal variances\n";
	std::cout<<"t-statistic:\t" + template_funcs::toString(t_stat) + "\n"; 
	std::cout<<"p-value:\t" +  template_funcs::toString(p_val) + "\n"; 
	if(p_val < 0.05){
		std::cout<<"Data have significantly different means\n\n";
	}
	else{
		std::cout<<"Data do not have significantly different means\n\n";
	}
}

void testing::perform_t_test_diff_mean_diff_var()
{
	// perform a Student's t-test on two distributions with different means and different variance
	// determine if the data // determine if the data have significantly different mean
	// R. Sheehan 11 - 9 - 2017

	std::string file1, file2; 
	std::vector<double> data1, data2;
	int n_data_1, n_data_2;
	double mean, var, t_stat, p_val; 

	bool loud = true; 

	//file1 = "Data_1.txt"; // mu = 0, sig = 1
	file1 = "Data_6.txt";  // mu = 1, sig = 1

	file2 = "Data_3.txt"; // mu = 1.5, sig = 0.5
	//file1 = "Data_4.txt";  // mu = 0, sig = 0.5
	//file2 = "Data_5.txt";  // mu = 0.01, sig = 0.5

	useful_funcs::read_into_vector(file1, data1, n_data_1, loud);

	statistic::avevar(data1, n_data_1, &mean, &var);

	std::cout<<"\nFile:\t" + file1 + "\n"; 
	std::cout<<"Size:\t" + template_funcs::toString(n_data_1) + "\n"; 
	std::cout<<"Average:\t" + template_funcs::toString(mean, 4) + "\n"; 
	std::cout<<"Variance:\t" + template_funcs::toString(var, 4) + "\n\n";

	useful_funcs::read_into_vector(file2, data2, n_data_2, loud);

	statistic::avevar(data2, n_data_2, &mean, &var);

	std::cout<<"\nFile:\t" + file2 + "\n"; 
	std::cout<<"Size:\t" + template_funcs::toString(n_data_2) + "\n"; 
	std::cout<<"Average:\t" + template_funcs::toString(mean, 4) + "\n"; 
	std::cout<<"Variance:\t" + template_funcs::toString(var, 4) + "\n\n";

	statistic::tutest(data1, n_data_1, data2, n_data_2, &t_stat, &p_val);

	std::cout<<"\nResult of Student's t-test, assuming unequal variances\n";
	std::cout<<"t-statistic:\t" + template_funcs::toString(t_stat) + "\n"; 
	std::cout<<"p-value:\t" +  template_funcs::toString(p_val) + "\n"; 
	if(p_val < 0.05){
		std::cout<<"Data have significantly different means\n\n";
	}
	else{
		std::cout<<"Data do not have significantly different means\n\n";
	}
}

void testing::perform_f_test_diff_mean_diff_var()
{
	// perform a Fischer F-test on two distributions with different means and different variance
	// determine if the data have significantly different variance
	// R. Sheehan 11 - 9 - 2017

	std::string file1, file2; 
	std::vector<double> data1, data2;
	int n_data_1, n_data_2;
	double mean, var, t_stat, p_val; 

	bool loud = true; 

	//file1 = "Data_1.txt"; // mu = 0, sig = 1
	//file2 = "Data_6.txt";  // mu = 1, sig = 1

	file1 = "Data_3.txt"; // mu = 1.5, sig = 0.5
	file2 = "Data_4.txt";  // mu = 0, sig = 0.5
	//file2 = "Data_5.txt";  // mu = 0.01, sig = 0.5

	useful_funcs::read_into_vector(file1, data1, n_data_1, loud);

	statistic::avevar(data1, n_data_1, &mean, &var);

	std::cout<<"\nFile:\t" + file1 + "\n"; 
	std::cout<<"Size:\t" + template_funcs::toString(n_data_1) + "\n"; 
	std::cout<<"Average:\t" + template_funcs::toString(mean, 4) + "\n"; 
	std::cout<<"Variance:\t" + template_funcs::toString(var, 4) + "\n\n";

	useful_funcs::read_into_vector(file2, data2, n_data_2, loud);

	statistic::avevar(data2, n_data_2, &mean, &var);

	std::cout<<"\nFile:\t" + file2 + "\n"; 
	std::cout<<"Size:\t" + template_funcs::toString(n_data_2) + "\n"; 
	std::cout<<"Average:\t" + template_funcs::toString(mean, 4) + "\n"; 
	std::cout<<"Variance:\t" + template_funcs::toString(var, 4) + "\n\n";

	statistic::ftest(data1, n_data_1, data2, n_data_2, &t_stat, &p_val);

	std::cout<<"\nResult of Fischer's f-test\n";
	std::cout<<"t-statistic:\t" + template_funcs::toString(t_stat) + "\n"; 
	std::cout<<"p-value:\t" +  template_funcs::toString(p_val) + "\n"; 
	if(p_val < 0.05){
		std::cout<<"Data have significantly different variances\n\n";
	}
	else{
		std::cout<<"Data do not have significantly different variances\n\n";
	}
}

void testing::perform_chsone_test()
{
	// perform chi-sq test to determine if the data in bins is from the expected distribution given by ebins
	// ebins = expected distribution
	// bins = actual data
	// chisq = chi-sq value
	// prob = probability that ebins and bins represent the same distribution
	// R. Sheehan 28 - 9 - 2017

	int NBINS = 10, NPTS = 2000; 

	int i, ibin, k; 
	long idum = -8; 
	double chisq, df, prob, x; 

	//df = degrees of freedom = NBINS

	// declare some vectors to hold the binned data
	std::vector<double> bins(NBINS,0.0); 
	std::vector<double> ebins(NBINS); 

	// generate exponentially distributed data bins
	// How is the width of each bin determined? 
	for(k=0; k<NPTS; k++){
		x = rng::expdev(&idum); // generate exponentially distributed random number
		ibin = static_cast<int>( x*NBINS/3.0 ); // determine bin into which x is to be placed
		//std::cout<<"x = "<<x<<", ibin = "<<ibin<<"\n"; 
		if(ibin < NBINS) ++bins[ibin]; // increment the number of counts in a bin
	}

	/*std::cout<<"The binned data is\n"; 
	for(i=0;i<NBINS;i++) std::cout<<bins[i]<<"\n"; 
	std::cout<<"\n";*/ 

	// generate binned data from a known distribution for comparison
	for(i=0; i<NBINS; i++){
		ebins[i] = 3.0*NPTS/NBINS*exp(-3.0*(i+0.5)/NBINS); 
	}

	// perform chi-sq test to determine if data if bins is from same distribution as data in ebins
	int knstrn = 0; 
	statistic::chsone(bins, ebins, NBINS, knstrn, &df, &chisq, &prob); 

	std::cout<<std::setw(15)<<"expected"<<std::setw(16)<<"observed"<<"\n"; 
	std::cout<<std::fixed<<std::setprecision(2); 
	for(i=0; i<NBINS; i++){
		std::cout<<std::setw(14)<<ebins[i]<<std::setw(16)<<bins[i]<<"\n"; 
	}
	std::cout<<"\n"<<std::setw(19)<<"chi-squared";
	std::cout<<std::setw(11)<<chisq<<"\n"; 
	std::cout<<std::setw(19)<<"probability:"<<std::setw(11)<<prob<<"\n"; 
	//std::cout<<std::setw(19)<<"df:"<<std::setw(11)<<df<<"\n"; 
}

void testing::perform_chstwo_test()
{
	// perform chi-sq test to determine if the data in bins1 is from the same distribution as by bins2
	// chisq = chi-sq value
	// prob = probability that ebins and bins represent the same distribution
	// R. Sheehan 28 - 9 - 2017

	int NBINS = 10, NPTS = 2000; 

	int i, ibin, k; 
	long idum = -17; 
	double chisq, df, prob, x; 

	//df = degrees of freedom = NBINS

	// declare some vectors to hold the binned data
	std::vector<double> bins1(NBINS,0.0); 
	std::vector<double> bins2(NBINS,0.0); 

	// generate exponentially distributed data bins
	// How is the width of each bin determined? 
	for(k=0; k<NPTS; k++){

		x = rng::expdev(&idum); // generate exponentially distributed random number
		ibin = static_cast<int>( x*NBINS/3.0 ); // determine bin into which x is to be placed
		if(ibin < NBINS) ++bins1[ibin]; // increment the number of counts in a bin

		x = rng::expdev(&idum); // generate exponentially distributed random number
		ibin = static_cast<int>( x*NBINS/3.0 ); // determine bin into which x is to be placed
		if(ibin < NBINS) ++bins2[ibin]; // increment the number of counts in a bin
	}

	/*std::cout<<"The binned data is\n"; 
	for(i=0;i<NBINS;i++) std::cout<<bins[i]<<"\n"; 
	std::cout<<"\n";*/ 

	// perform chi-sq test to determine if data from bins1 is from same distribution as data in bins2
	int knstrn = 0; 
	statistic::chstwo(bins1, bins2, NBINS, knstrn, &df, &chisq, &prob); 

	std::cout<<std::setw(15)<<"data set 1"<<std::setw(16)<<"data set 2"<<"\n"; 
	std::cout<<std::fixed<<std::setprecision(2); 
	for(i=0; i<NBINS; i++){
		std::cout<<std::setw(14)<<bins1[i]<<std::setw(16)<<bins2[i]<<"\n"; 
	}
	std::cout<<"\n"<<std::setw(19)<<"chi-squared";
	std::cout<<std::setw(11)<<chisq<<"\n"; 
	std::cout<<std::setw(19)<<"probability:"<<std::setw(11)<<prob<<"\n"; 
	//std::cout<<std::setw(19)<<"df:"<<std::setw(11)<<df<<"\n"; 
}

double testing::ksfunc(double x)
{
	//return 1.0 - probability::erffc(x/sqrt(2.0)); 
	return probability::erff(x/sqrt(2.0)); 
}

void testing::perform_ksone_test()
{
	// routine for examining the Kolmorgoriv-Smirnov test
	// ksone use K-S criterion to compare single data set to an expected distribution
	// kstwo use K-S criterion to compare two data sets
	// this program creates data sets with Gaussian distributions and stepwise increasing variance
	// and compares their cumulative distribution function to the expected result for Gaussian 
	// distribution of unit variance, i.e the error function erff
	// Increasing variance in the test distribution should reduce the likelihood that it was drawn from
	// the same distribution represented by the comparison function
	// R. Sheehan 28 - 9 - 2017

	int NPTS = 1000; 
	int i, j; 
	long idum = (-5); 
	double epsilon = 0.1, d, factr, prob, varnce;

	std::vector<double> data(NPTS); 

	std::cout<<std::setw(19)<<"variance ratio"<<std::setw(17)<<"k-s statistic";
	std::cout<<std::setw(16)<<"probability"<<"\n\n"; 
	std::cout<<std::fixed<<std::setprecision(7);
	for(i=0; i<11; i++){
		varnce = 1.0+i*epsilon;
		factr = sqrt(varnce); 
		for(j=0; j<NPTS; j++){
			data[j] = factr*fabs(rng::gasdev(&idum)); 
		}
		statistic::ksone(data, NPTS, ksfunc, &d, &prob);
		std::cout<<std::setw(16)<<varnce<<std::setw(17)<<d;
		std::cout<<std::setw(17)<<prob<<"\n"; 
	}
}

void testing::perform_kstwo_test()
{
	// routine for examining the Kolmorgoriv-Smirnov test
	// ksone use K-S criterion to compare single data set to an expected distribution
	// kstwo use K-S criterion to compare two data sets
	// kstwo compares the cumulative distribution functions of two unbinned data sets data1 and data2
	// In this program both data sets are Gaussian distributions but data2 is given a stepwise increase in variance

	int N1 = 1000, N2 = 2000, i, j;
	long idum = (-1357); 
	double epsilon = 0.1, d, factr, prob, varnce; 
	std::vector<double> data1(N1), data2(N2); 

	for(i=0; i<N1; i++) data1[i] = rng::gasdev(&idum); 

	std::cout<<std::setw(18)<<"vaariance ratio"<<std::setw(16)<<"k-s statistic";
	std::cout<<std::setw(15)<<"probability\n\n"; 
	std::cout<<std::fixed<<std::setprecision(6);
	for(i=0; i<11; i++){
		varnce = 1.0+i*epsilon; 
		factr = sqrt(varnce); 
		for(j=0; j<N2; j++) data2[j] = factr*rng::gasdev(&idum); 
		statistic::kstwo(data1, N1, data2, N2, &d, &prob); 
		std::cout<<std::setw(15)<<varnce<<std::setw(16)<<d;
		std::cout<<std::setw(16)<<prob<<"\n"; 
	}

}

void testing::read_tab1_data(std::string &header, std::vector<std::string> &months, std::vector<std::string> &accident_types, std::vector<std::vector<int>> &numbers)
{
	// Read the data stored in the file TABLE1.DAT
	// store the names of the months in the vector months
	// store the types of accidents listed in vector accident_types
	// store the numbers of accidents in the vector numbers, columns = month, row = accident type
	// R. Sheehan 11 - 1 - 2018

	try{
		std::string filename = "TABLE1.DAT";

		std::ifstream the_file; 
		the_file.open(filename, std::ios_base::in); 
		if(the_file.is_open()){

			std::string line, item; 

			int line_num = 0; 
			while( std::getline(the_file, line) ){

				//extract the header from the file
				if(line_num==0)header = line + '\n'; 
				if(line_num==1)header += line; 

				// parse the stream to extract the months
				if(line_num == 2){					
					std::istringstream linestream(line);
					while(std::getline(linestream, item, ' ')){
						if( !(item.find("Month:") != std::string::npos) ) { // only save the string if it is the name of a month
							std::remove_if(item.begin(), item.end(), isspace); // remove all white space from the string
							if( !item.empty() ) months.push_back(item); 
						}
					}
				}
				
				// parse the stream to extract types and numbers of accident for each month
				if(line_num>3){
					std::istringstream linestream(line);
					int item_cnt = 0; 
					std::vector<int> acc_per_type; // vector to hold the numbers of accidents of each type
					while(std::getline(linestream, item, ' ')){					
						std::remove_if(item.begin(), item.end(), isspace); // remove all white space from the string
						if(item_cnt == 0 && !item.empty() ) accident_types.push_back( item ); // save the names of the types of accidents
						if(item_cnt > 0 && !item.empty() ) acc_per_type.push_back( std::stoi(item) ); // convert the number from string to int, atof(item.c_str()) will also work
						item_cnt++; 
					}
					numbers.push_back(acc_per_type); 
					acc_per_type.clear(); 
				}

				line_num++; 
			}

			the_file.close(); 
		}
		else{
			std::string reason = "Error: Test_Routines::read_tab1_data\nCannot open file: " + filename + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
	}
}

void testing::read_tab2_data(std::string &header, std::vector<std::string> &months, std::vector<std::string> &cities, std::vector<std::string> &states, std::vector<double> &state_averages, std::vector<double> &state_latitudes, std::vector<std::vector<double>> &numbers)
{
	// Read the data stored in the file TABLE2.DAT
	// store the names of the months in the vector months
	// store the names of the states in the vector states
	// store the state averages in the vector state_averages
	// store the state latitudes in the vector state_latitude
	// store the average monthly solar radiation in the vector numbers, columns = month, row = state
	// R. Sheehan 29 - 3 - 2018

	try{
		std::string filename = "TABLE2.DAT";

		std::ifstream the_file; 
		the_file.open(filename, std::ios_base::in); 
		if(the_file.is_open()){
			std::string line, item; 

			int line_num = 0; 
			while( std::getline(the_file, line) ){

				//extract the header from the file
				if(line_num==0)header = line + '\n'; 
				if(line_num==1)header += line; 

				// parse the stream to extract the months
				if(line_num == 2){					
					std::istringstream linestream(line);
					while(std::getline(linestream, item, ' ')){
						bool c1 = !(item.find("Month:") != std::string::npos); 
						bool c2 = !(item.find("ave") != std::string::npos);
						bool c3 = !(item.find("lat") != std::string::npos);
						if( c1 && c2 && c3 ) { // only save the string if it is the name of a month
							std::remove_if(item.begin(), item.end(), isspace); // remove all white space from the string
							if( !item.empty() ) months.push_back(item); 
						}
					}
				}

				// parse the stream to extract state, numbers for each month, state averages and state latitudes
				if(line_num>3){
					std::istringstream linestream(line);
					int item_cnt = 0, mth_cnt = 0, num_mths = 12, end_rad_per_state, avg_pos; 
					bool get_avg, get_lat; 
					std::vector<double> rad_per_state; // vector to hold the numbers of accidents of each type
					while(std::getline(linestream, item, ' ')){							
						std::remove_if(item.begin(), item.end(), isspace); // remove all white space from the string
						if(item_cnt == 0 && !item.empty() ){
							cities.push_back( item ); // save the names of the cities
							end_rad_per_state = -1; avg_pos = -1;
							get_avg = false; get_lat = false; 
						}
						if(item_cnt == 1 && !item.empty() ) states.push_back( item ); // save the names of the cities						
						if(item_cnt > 1 && !item.empty() && (int)(rad_per_state.size())< num_mths){
							rad_per_state.push_back( std::stof(item) ); // convert the number from string to int, atof(item.c_str()) will also work
							mth_cnt++; 
						}
						if(mth_cnt == 12){
							end_rad_per_state = item_cnt; // mark position where all month data per line is stored in vector
							mth_cnt = 0; 
							get_avg = true; get_lat = false; 
						}
						// extract state averages and latitudes						
						if(get_avg && item_cnt > end_rad_per_state && !item.empty()){
							state_averages.push_back(std::stof(item));
							avg_pos = item_cnt;
							get_avg = false; get_lat = true; 
						}
						// extract state latitudes
						if(get_lat && item_cnt > avg_pos && !item.empty()){
							state_latitudes.push_back(std::stof(item));
							get_avg = false; get_lat = false; 
						}
						item_cnt++; 
					}
					numbers.push_back(rad_per_state); 
					rad_per_state.clear(); 
				}

				line_num++; 
			}

			the_file.close(); 
		}
		else{
			std::string reason = "Error: Test_Routines::read_tab2_data\nCannot open file: " + filename + "\n";
			throw std::invalid_argument(reason);
		}

	
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
	}
}

void testing::perform_con_tab_test()
{
	// Test the routines for computing the measures of association for a contingency table

	// Read the data from the file
	std::string Head; 

	std::vector<std::string> Mths; 

	std::vector<std::string> Acc; 

	std::vector< std::vector<int> > Num; 
	
	testing::read_tab1_data(Head, Mths, Acc, Num);

	// perform test based on Cramer's V and CC
	int nrows, ncols; 
	double chisq, df, prob, cramrv, ccc; 

	nrows = static_cast<int>(Num.size()); 
	ncols = static_cast<int>(Num[0].size()); 

	statistic::cntab1(Num, nrows, ncols, &chisq, &df, &prob, &cramrv, &ccc);

	std::cout<<"Results from Contingency Table Analysis\n";
	std::cout<<std::fixed<<std::setprecision(4);
	std::cout<<std::endl<<std::setw(20)<<"chi-squared";
	std::cout<<std::setw(20)<<chisq<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"degrees of freedom";
	std::cout<<std::setw(20)<<df<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"probability";
	std::cout<<std::setw(20)<<prob<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"cramer-V";
	std::cout<<std::setw(20)<<cramrv<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"contingency coeff.";
	std::cout<<std::setw(20)<<ccc<<std::endl<<std::endl;

	// perform test based on measures of entropy
	double h, hx, hy, hygx, hxgy, uygx, uxgy, uxy; 

	statistic::cntab2(Num, nrows, ncols, &h, &hx, &hy, &hygx, &hxgy, &uygx, &uxgy, &uxy); 

	std::cout<<"Results from entropy-based Contingency Table Analysis\n";
	std::cout<<std::fixed<<std::setprecision(4);
	std::cout<<std::endl<<std::setw(20)<<"entropy of table";
	std::cout<<std::setw(20)<<h<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"entropy of x";
	std::cout<<std::setw(20)<<hx<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"entropy of y";
	std::cout<<std::setw(20)<<hy<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"entropy of y|x";
	std::cout<<std::setw(20)<<hygx<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"entropy of x|y";
	std::cout<<std::setw(20)<<hxgy<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"dependency of y|x";
	std::cout<<std::setw(20)<<uygx<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"dependency of x|y";
	std::cout<<std::setw(20)<<uxgy<<std::endl;
	std::cout<<std::endl<<std::setw(20)<<"symmetrical depend.";
	std::cout<<std::setw(20)<<uxy<<std::endl;
}

void testing::perform_pearson_test()
{
	// routine for testing the calculation of Pearson correlation coefficient
	// this program examines two variables to find linear correlations
	// It returns linear correlation coefficient r, a probability of correlation prob and Fisher's z
	// This sample program sets up data dose and spore which show hypothetical spore count from plants 
	// exposed to various levels of gamma rays
	// R. Sheehan 29 - 3 - 2018

	int N; 
	double prob, r, z; 

	std::string dose_file = "Dose.txt";
	std::string spore_file = "Spore.txt";

	std::vector<double> dose; 
	std::vector<double> spore; 

	// read the data from the files
	useful_funcs::read_into_vector(dose_file, dose, N); 
	useful_funcs::read_into_vector(spore_file, spore, N); 

	// print the data to screen
	std::cout<<"Effect of Gamma Rays on Marigolds\n"; 
	std::cout<<std::setw(16)<<"Count Rate (cpm)"; 
	std::cout<<std::setw(24)<<"Pollen Index\n";
	std::cout<<std::fixed<<std::setprecision(2);
	for(int i = 0; i < N; i++){
		std::cout<<std::setw(10)<<dose[i]<<std::setw(26)<<spore[i]<<"\n"; 
	}

	// Compute Pearson's r correlation coefficient and related probabilities
	statistic::pearsn(dose, spore, N, &r, &prob, &z); 

	std::cout<<"\n"<<std::setw(30)<<"PEARSN"<<std::setw(17)<<"Expected\n"; 
	std::cout<<"Corr. Coeff."<<std::setw(9)<<" ";
	std::cout<<std::fixed<<std::setprecision(6);
	std::cout<<std::setw(10)<<r<<std::setw(16)<<0.9069586<<"\n";
	std::cout<<"Probability"<<std::setw(10)<<" ";
	std::cout<<std::setw(10)<<prob<<std::setw(16)<<0.2926505e-3<<"\n";
	std::cout<<"Fisher's z"<<std::setw(11)<<" ";
	std::cout<<std::setw(10)<<z<<std::setw(16)<<1.510110<<"\n";

	// compute the significance level at which r differs from r_{true}, assuming z^{bar}=0
	double test_val = fabs(z)*sqrt(0.5*((double)(N)-3.0)); 
	double alt_p_value = probability::erffc(test_val); 

	std::cout<<"Significance Level at which r differs from r-true: "<<alt_p_value<<"\n"; 
}

void testing::perform_spearman_test()
{
	// Rank order correlation is done with spear to compare two distributions data1 and data2 for correlation
	// Correlations are reported in terms of d and rs
	// d = sum squared difference in ranks
	// rs = Spearman's rank order correlation parameter
	// This program looks at the solar radiation data in TABLE2.DAT and looks for correlations between columns of the table
	// with each column considered a separate data set
	// R. Sheehan 29 - 3 - 2018

	// containers for storing the data
	int i, j; 
	int NDAT = 20, NMON = 12; 
	double d, probd, probrs, rs, zd; 
	std::string Head; 
	std::vector<std::string> Mths; 	
	std::vector<std::string> Cities;
	std::vector<std::string> States;
	std::vector<double> State_Averages;
	std::vector<double> State_Latitudes;
	std::vector< std::vector<double> > Rays;
	std::vector<double> data1(NDAT); 
	std::vector<double> data2(NDAT); 

	// read the data from the file
	read_tab2_data(Head, Mths, Cities, States, State_Averages, State_Latitudes, Rays);

	// Check temperature correlations between different months
	std::cout<<"Are sunny summer places also sunny winter places?\n";
	std::cout<<"Check correlation of samples US solar radiation\n";
	std::cout<<"Compare July with other months\n\n"; 
	std::cout<<"month"<<std::setw(10)<<"d"<<std::setw(15)<<"std. dev.";
	std::cout<<std::setw(12)<<"probd"<<std::setw(16)<<"spearman-r";
	std::cout<<std::setw(11)<<"probrs\n\n";
	
	// store july data in data1
	for(i=0;i<NDAT;i++){
		data1[i] = Rays[i][0]; 
	}

	// loop over correlation calculation for remaining months
	for(j=0; j<NMON; j++){		
		
		// store data for month j in data2
		for(i=0; i<NDAT; i++){
			data2[i] = Rays[i][j]; 
		}

		// perform spear rank order correlation test
		statistic::spear(data1, data2, NDAT, &d, &zd, &probd, &rs, &probrs); 

		// output computed results
		std::cout<<std::fixed<<std::setprecision(2); 
		std::cout<<std::setw(4)<<Mths[j]<<std::setw(13)<<d; 
		std::cout<<std::fixed<<std::setprecision(6); 
		std::cout<<std::setw(13)<<zd<<std::setw(13)<<probd<<std::setw(14)<<rs; 
		std::cout<<std::setw(13)<<probrs<<"\n"; 
	}
}

void testing::rank_order_test()
{
	// This program looks at the solar radiation data in TABLE2.DAT and computes the rank order for the data therein
	// Each column of the solar flux data is replaced by the rank order of its entries
	// R. Sheehan 12 - 4 - 2018

	// Containers for holding the data
	int NDAT = 20, NMON = 12;
	int i, j; 
	std::string Head;
	std::vector<std::string> Mths;
	std::vector<std::string> Cities;
	std::vector<std::string> States;
	std::vector<double> State_Averages;
	std::vector<double> State_Latitudes;
	std::vector< std::vector<double> > Rays;
	std::vector<double> data(NDAT);
	std::vector<double> order(NDAT); 
	std::vector<double> s(NMON); 

	// read the data from the file
	read_tab2_data(Head, Mths, Cities, States, State_Averages, State_Latitudes, Rays);

	// replace the solar flux in each column by rank order
	for (j = 0; j < NMON; j++) {
		for (i = 0; i < NDAT; i++) {
			data[i] = Rays[i][j];
			order[i] = i;
		}
		template_funcs::sort2(data, order);
		statistic::crank(NDAT, data, &s[j]);
		for (i = 0; i < NDAT; i++) {
			Rays[(int)(0.5 + order[i])][j] = data[i];
		}
	}
	std::cout << Head << "\n\n" << std::setw(16) <<" ";
	for (i = 0; i < NMON; i++) {
		std::cout << std::setw(5) << Mths[i]; 
	}
	std::cout << "\n\n";
	for (i = 0; i < NDAT; i++) {
		std::cout << Cities[i];
		for (j = 0; j < NMON; j++) {
			std::cout << std::setw(5) << (int)(0.5 + Rays[i][j]);
		}
		std::cout << "\n"; 
	}
}

void testing::perform_kendall_test_1()
{
	// This program looks for pair correlations in five random number generators
	// The data are of the same rank
	// it tests for randomness by seeing if two consecutive numbers from the generator have a monotonic correlation
	// RNG generate 200 pairs of random numbers, then kendl1 tests for correlation of the pairs
	// A chart is made showing Kendall's tau, std. dev from null hypothesis and propability
	// To improve the test increase NDAT and vary the seed value idum
	// R. Sheehan 12 - 4 - 2018

	int NDAT = 500; 
	int i, j; 
	double prob, tau, z; 

	std::vector<std::string> txt = {"RAN0", "RAN1", "RAN2"}; //{"RAN0", "RAN1", "RAN2", "RAN3", "RAN4"};
	std::vector<double> data1(NDAT); 
	std::vector<double> data2(NDAT); 

	// Look for correlations in "RAN0", "RAN1", "RAN2"
	std::cout << "\nPair correlations of RAN0, RAN1, RAN2\n\n"; 
	std::cout << std::setw(9) << "\nProgram" << std::setw(18) << "Kendall tau"; 
	std::cout << std::setw(17) << "Deviation" << std::setw(19) << "Probability\n";
	std::cout << std::fixed << std::setprecision(6); 
	for (i = 0; i < 3; i++) {
		long idum = (-1357); 
		// generate pairs of random numbers using the different rng
		for (j = 0; j < NDAT; j++) {
			if (i == 0) {
				data1[j] = rng::ran0(&idum);
				data2[j] = rng::ran0(&idum);
			}
			else if(i == 1) {
				data1[j] = rng::ran1(&idum);
				data2[j] = rng::ran1(&idum);
			}
			else if (i == 2) {
				data1[j] = rng::ran2(&idum);
				data2[j] = rng::ran2(&idum);
			}
			else {
				data1[j] = rng::ran2(&idum);
				data2[j] = rng::ran2(&idum);
			}
		}
		statistic::kendl1(data1, data2, NDAT, &tau, &z, &prob); 
		std::cout << std::setw(8) << txt[i] << std::setw(18) << tau;
		std::cout << std::setw(18) << z << std::setw(18) << prob<<"\n";
	}
}

void testing::perform_kendall_test_2()
{
	// Program perform Kendall tau test on a contingency table of data of random bit sequences
	// Contingency table is prepared inside the program
	// Program checks the sequences by breakig them into groups of three bits each. Each group is treated as a 3-bit binary number
	// Two consecutive groups then act as indices in an 8*8 contingency table that records how many times each possible sequence occurs
	// For each random bit generator NDAT samples are taken. Kendall tau and associated probability are computed using kendl2
	// Program is testing whether or not the 3-bit numbers tend to be followed by others of their own kind
	// Roughly this program tests whether groups of ones or zeroes come in more groups than they should
	// R. Sheehan 12 - 4 -2018

	int NDAT = 1000, IP = 8, JP = 8; 
	int ifunc, k, l, m, n, twoton; 
	unsigned long iseed; 
	double prob, tau, z; 

	std::vector<std::string> txt = { "000", "001", "010", "011", "100", "101", "110", "111" };
	std::vector< std::vector< double > > tab(IP); 
	for (k = 0; k < IP; k++) {
		std::vector<double> x(JP); 
		tab[k] = x; 
	}

	// look for ones after zeroes in IRBIT1 and IRBIT2 sequences
	std::cout << "Are ones followed by zeroes and vice-versa?\n"; 
	std::cout << std::fixed << std::setprecision(6); 
	for (ifunc = 0; ifunc < 2; ifunc++) {
		iseed = 2468; 
		
		if (ifunc == 0) {
			std::cout << "test of irbit1:\n"; 
		}
		else {
			std::cout << "test of irbit2:\n"; 
		}
		
		// zero the array holding the data
		for (k = 0; k < IP; k++) {
			for (l = 0; l < JP; l++) {
				tab[k][l] = 0.0; 
			}
		}
		
		// Fill the array with data
		for (m = 1; m <= NDAT; m++) {
			k = 0; 
			twoton = 1; 
			for (n = 0; n < 3; n++) {
				if (ifunc == 0) {
					k += rng::irbit1(&iseed)*twoton; 
				}
				else {
					k += rng::irbit2(&iseed)*twoton; 
				}
				twoton *= 2; 
			}
			l = 0; 
			twoton = 1;
			for (n = 0; n < 3; n++) {
				if (ifunc == 0) {
					l += rng::irbit1(&iseed)*twoton;
				}
				else {
					l += rng::irbit2(&iseed)*twoton;
				}
				twoton *= 2;
			}
			++tab[k][l]; 
		}
		
		// Compute Kentall tau statistic
		statistic::kendl2(tab, IP, JP, &tau, &z, &prob);
		
		// print data to screen
		std::cout << "    ";
		for (n = 0; n < 8; n++) {
			std::cout << std::setw(6) << txt[n];
		}
		std::cout << "\n"; 
		for (n = 0; n < 8; n++) {
			std::cout << std::setw(3) << txt[n]; 
			for (m = 0; m < 8; m++)
				std::cout << std::setw(6) << (int)(0.5 + tab[n][m]);
			std::cout << "\n"; 
		}
		std::cout << "\n" << std::setw(17) << "Kendall tau:";
		std::cout << std::setw(15) << "deviation" << std::setw(17) << "probability\n";
		std::cout << std::setw(15) << tau << std::setw(16) << z; 
		std::cout << std::setw(16) << prob << "\n\n";
	}
}

void testing::lin_fit_test()
{
	// functions makes a linear fit to a set of data points with noise
	// underlying model is linear of the form y = a + b x
	// where a = +1 and b = -2, noise is assumed to be in y and not in x
	// R. Sheehan 22 - 6 - 2018

	int NPT = 100; 
	double SPREAD = 1.0; 
	int i, mwt = 0; 
	long idum = (-117); 
	double a, b, chi2, q, siga, sigb; 
	std::vector<double> x(NPT, 0.0); 
	std::vector<double> y(NPT, 0.0);
	std::vector<double> sig(NPT, 0.0);

	// fill the vectors with data
	for (i = 0; i < NPT; i++) {
		x[i] = 0.1*(i + 1); 
		y[i] = -2.0*x[i] + 1.0 + SPREAD*rng::gasdev(&idum); 
		sig[i] = SPREAD; 
	}

	// Compute the linear fits with / without assumed spread on y values
	std::cout << std::fixed << std::setprecision(6); 
	for (i = 0; i < 2; i++) {
		// compute the linear fit parameters
		fit::lin_fit(x, y, NPT, sig, mwt, &a, &b, &siga, &sigb, &chi2, &q); 

		if (mwt)
			std::cout << "\nIncluding standard deviations\n"; 
		else
			std::cout << "\nIgnoring standard deviations\n";
		std::cout << std::setw(12) << "a = " << std::setw(10) << a; 
		std::cout << std::setw(19) << "uncertainty: " << std::setw(10) << siga << "\n";
		std::cout << std::setw(12) << "b = " << std::setw(10) << b;
		std::cout << std::setw(19) << "uncertainty: " << std::setw(10) << sigb << "\n";
		std::cout << std::setw(19) << "chi-squared = " << std::setw(15) << chi2 << "\n";
		std::cout << std::setw(19) << "nu = " << std::setw(15) << NPT - 2 << "\n";
		std::cout << std::setw(19) << "chi-squared / nu = " << std::setw(15) << chi2 / (NPT - 2) << "\n";
		std::cout << std::setw(23) << "goodness-of-fit = " << std::setw(11) << q << "\n";
		mwt++; 
	}
	std::cout << "\n"; 
}

void testing::gaussj_test()
{
	// Ensure the Gauss elimination scheme is working correctly
	// R. Sheehan 22 - 6 - 2018

	int nrows = 4, ncols = 4, ncolsb = 1; 
	std::vector<std::vector<double>> A; 
	std::vector<std::vector<double>> b;

	A = vecut::array_2D(nrows, ncols); 
	b = vecut::array_2D(nrows, ncolsb);

	// Fill the matrix A
	A[0][0] = 1; A[0][1] = -1; A[0][2] = 2; A[0][3] = -1; 
	A[1][0] = 2; A[1][1] = -2; A[1][2] = 3; A[1][3] = -3;
	A[2][0] = 1; A[2][1] =  1; A[2][2] = 1; A[2][3] = 0;
	A[3][0] = 1; A[3][1] = -1; A[3][2] = 4; A[3][3] = 3;

	// Fill the vector B
	b[0][0] = -8; 
	b[1][0] = -20; 
	b[2][0] = -2; 
	b[3][0] = 4;

	std::cout << "Matrix A is\n";
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++)
			std::cout << A[i][j] << " ";
		std::cout << "\n";
	}
	std::cout << "\nVector b is\n";
	for (int i = 0; i < nrows; i++) {
		std::cout << b[i][0] << "\n";
	}

	lin_alg::gaussj(A, nrows, b, ncolsb); 

	std::cout << "\nInverse of A is\n"; 
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < ncols; j++)
			std::cout << A[i][j] << " "; 
		std::cout << "\n"; 
	}
	std::cout << "\nSolution of system of equations is\n";
	for (int i = 0; i < nrows; i++) {
		std::cout << b[i][0] << "\n";
	}

	A.clear(); b.clear(); 
}

void testing::funcs(double &x, std::vector<double> &afunc, int &ma)
{
	// function for producing sample basis function which can be used to fit data
	// 
	try {
		if (ma > 2 && ma == (int)(afunc.size())) {
			int i;
			afunc[0] = 1.0;
			afunc[1] = x;
			for (i = 2; i < ma; i++) afunc[i] = sin((i + 1)*x);
		}
		else{
			std::string reason = "Error: testing::funcs()\n"; 
			reason += "ma = " + template_funcs::toString(ma) + ", afuncs.size() == " + template_funcs::toString(afunc.size()) + "\n"; 
			throw std::runtime_error(reason); 
		}
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what(); 
	}
}

void testing::fpoly(double &x, std::vector<double> &afunc, int &ma)
{
	// General polynomial with known coefficients	
	try {		
		if (ma > 2 && ma == (int)(afunc.size())) {
			int j;
			afunc[0] = 1.0;
			for (j = 1; j < ma; j++) afunc[j] = afunc[j - 1] * x;
		}
		else {
			std::string reason = "Error: testing::fpoly()\n";
			reason += "ma = " + template_funcs::toString(ma) + ", afuncs.size() == " + template_funcs::toString(afunc.size()) + "\n";
			throw std::runtime_error(reason);
		}
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
	}
}

void testing::fleg(double &x, std::vector<double> &afunc, int &ma)
{
	// Compute coefficients of legendre polynomial using recurrence relation

	try {
		if (ma > 2 && ma == (int)(afunc.size())) {
			int j;
			double twox, f2, f1, d; 
			afunc[0] = 1.0;
			afunc[1] = x; 
			twox = 2.0*x; 
			f2 = x; 
			d = 1.0; 
			for (j = 2; j < ma; j++) { 
				f1 = d++; 
				f2 += twox; 
				afunc[j] = (f2*afunc[j-1]-f1*afunc[j-2]) / d;
			}
		}
		else {
			std::string reason = "Error: testing::fleg()\n";
			reason += "ma = " + template_funcs::toString(ma) + ", afuncs.size() == " + template_funcs::toString(afunc.size()) + "\n";
			throw std::runtime_error(reason);
		}
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what();
	}
}

void testing::lst_sqr_test()
{
	// Make a polynomial fit to a set of data
	// 
	int i, j, NPT = 100, NTERM = 5, nu = NPT - NTERM; 
	double SPREAD = 0.1, chisq, q; 
	long idum = (-911); 

	std::vector<int> ia(NTERM); 
	std::vector<double> a(NTERM), x(NPT), y(NPT), sig(NPT); 
	std::vector<std::vector<double>> covar = vecut::array_2D(NTERM, NTERM);

	// Generate data for the fit
	for (i = 0; i < NPT; i++) {
		x[i] = 0.1*(i + 1); 
		funcs(x[i], a, NTERM); 
		y[i] = 0.0; 
		for (j = 0; j < NTERM; j++) y[i] += (j + 1)*a[j]; 
		y[i] += SPREAD * rng::gasdev(&idum); 
		sig[i] = SPREAD; 
	}

	// store which parameters are to be fitted
	for (i = 0; i < NTERM; i++) ia[i] = 1; 

	// perform least squares fit
	fit::lfit(x, y, sig, NPT, a, ia, NTERM, covar, &chisq, funcs); 

	q = probability::gammq(0.5*nu, 0.5*(chisq)); // goodness of fit

	// output results
	std::cout << "\n" << std::setw(11) << "parameter"; 
	std::cout << std::setw(22) << "uncertainty\n"; 
	std::cout << std::fixed << std::setprecision(6); 
	for (i = 0; i < NTERM; i++) {
		std::cout << "a[" << i << "] = " << std::setw(8) << a[i]; 
		std::cout << std::setw(13) << sqrt(covar[i][i]) << "\n"; 
	}
	std::cout << "chi-squared = " << std::setw(12) << chisq << "\n";
	std::cout << "nu = " << nu << "\n"; 
	std::cout << "chi_squared / nu = " << chisq / nu << "\n";
	std::cout << "goodness-of-fit = " << q << "\n\n"; 
	std::cout << "full covariance matrix\n";
	std::cout << std::scientific << std::setprecision(4); 
	for (i = 0; i < NTERM; i++) {
		for (j = 0; j < NTERM; j++)
			std::cout << std::setw(15) << covar[i][j]; 
		std::cout << "\n"; 
	}

	std::cout << "\nCheck the results of restricting the fit parameters...\n"; 

	// Check the results of restricting the fit parameters

	// store which parameters are to be fitted
	for (i = 1; i < NTERM; i+=2) ia[i] = 0;

	// perform least squares fit
	fit::lfit(x, y, sig, NPT, a, ia, NTERM, covar, &chisq, funcs);

	nu = NPT - 3;
	q = probability::gammq(0.5*nu, 0.5*(chisq)); // goodness of fit

	// output results
	std::cout << "\n" << std::setw(11) << "parameter";
	std::cout << std::setw(22) << "uncertainty\n";
	std::cout << std::fixed << std::setprecision(6);
	for (i = 0; i < NTERM; i++) {
		std::cout << "a[" << i << "] = " << std::setw(8) << a[i];
		std::cout << std::setw(13) << sqrt(covar[i][i]) << "\n";
	}
	std::cout << "chi-squared = " << std::setw(12) << chisq << "\n";
	std::cout << "nu = " << nu << "\n";
	std::cout << "chi_squared / nu = " << chisq / nu << "\n";
	std::cout << "goodness-of-fit = " << q << "\n\n";
	std::cout << "full covariance matrix\n";
	std::cout << std::scientific << std::setprecision(4);
	for (i = 0; i < NTERM; i++) {
		for (j = 0; j < NTERM; j++)
			std::cout << std::setw(15) << covar[i][j];
		std::cout << "\n";
	}
	std::cout << "\n";
}

void testing::fgauss(double x, std::vector<double> &a, double *y, std::vector<double> &dyda, int &na)
{
	// y(x;a) is the sum of na/3 Gaussians. The amplitude, centre and width of the Gaussians are stored in a
	// a[i] = amp_{k}, a[i+1] = centre_{k}, a[i+2] = width_{k}
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// k = 1..na/3
	// R. Sheehan 12 - 7 - 2018

	try {
		int i;
		double fac, ex, arg;

		*y = 0.0;
		for (i = 0; i < na - 1; i += 3) {
			arg = (x - a[i + 1]) / a[i + 2];
			ex = exp(-arg * arg);
			fac = a[i] * ex*2.0*arg;
			*y += a[i] * ex;
			dyda[i] = ex;
			dyda[i + 1] = fac / a[i + 2];
			dyda[i + 2] = fac * arg / a[i + 2];
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}
}

double testing::convert_dBm_to_mW(double dBm_val)
{
	// convert a dBm power reading to a mW power reading
	return pow(10.0, (dBm_val / 10.0));
}

double testing::convert_mW_to_dBm(double mW_val)
{
	// convert a mW power reading to dBm power reading
	if (mW_val > 1.0e-10) {
		return 10.0 * log10(mW_val);
	}
	else {
		return -100.0;
	}
}

void testing::non_lin_fit_test()
{
	// routine which determines if the Levenberg-Margquardt method has been implemented correctly
	// R. Sheehan 13 - 7 - 2018

	// Going to use the multi-peak Gaussian with noise as the test program

	// Can you update this code to provide an estimate of the error associated with each fit parameter? 
	// R. Sheehan 14 - 11 - 2019

	// Generate data to use in the fit process
	int npts, npars; 
	long idum = (-1011); 
	double spread, xlow, xhigh, deltax, xpos, yval;

	xlow = 0.0; xhigh = 10.0; deltax = 0.02; 
	npts = (int)( 1.0 + ( (xhigh - xlow) / deltax) ); 

	std::vector<double> xdata(npts, 0.0); 
	std::vector<double> ydata(npts, 0.0);
	std::vector<double> sigdata(npts, 0.0);

	npars = 6; 
	std::vector<double> a(npars, 0.0); 
	
	std::vector<double> dyda(npars, 0.0);

	// data stored in array a in the form a[i] = amp_{k}, a[i+1] = centre_{k}, a[i+2] = width_{k}
	a[0] = 5.0; a[1] = 3.0; a[2] = 2.0; // Gaussian 1
	a[3] = 3.0; a[4] = 7.0; a[5] = 1.0; // Gaussian 2

	spread = 0.01; 
	xpos = xlow; 
	for (int i = 0; i < npts; i++) {

		fgauss(xpos, a, &yval, dyda, npars); // evaluate the Gaussian function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		sigdata[i] = spread * yval;

		xpos += deltax;
	}

	// Perform the best it search for the data set
	int ITMAX = 10; 

	double TOL = 0.001; 
	double chisq = 0.0; 

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars); 
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	// Good guesses
	a_guess[0] = 4.5; a_guess[1] = 2.2; a_guess[2] = 2.8;
	a_guess[3] = 2.5; a_guess[4] = 6.3; a_guess[5] = 1.4;

	// Bad guesses
	/*a_guess[0] = 4.0; a_guess[1] = -2.2; a_guess[2] = 1.8;
	a_guess[3] = 0.5; a_guess[4] = -3.3; a_guess[5] = 0.4;*/

	// Can you add something that tells the user that the fit is good or bad?
	// How to interpret the chisq value correctly in terms of goodness of fit?

	// Find three parameters
	//ia[0] = 0; ia[1] = 0; ia[2] = 0; 
	//a_guess[0] = 5.0; a_guess[1] = 2.0; a_guess[2] = 3.0; // these parameters are fixed
	//ia[3] = 1; ia[4] = 1; ia[5] = 1;
	//a_guess[3] = 2.5; a_guess[4] = 4.3; a_guess[5] = 1.4; // these are the guesses for the parameters being sought

	// Find four parameters
	//ia[0] = 0; ia[1] = 1; ia[2] = 1;
	//a_guess[0] = 5.0; a_guess[1] = 2.4; a_guess[2] = 1.2; // these parameters are fixed
	//ia[3] = 1; ia[4] = 0; ia[5] = 1;
	//a_guess[3] = 2.5; a_guess[4] = 6.0; a_guess[5] = 1.4; // these are the guesses for the parameters being sought

	// Find two parameters
	/*ia[0] = 1; ia[1] = 0; ia[2] = 0;
	a_guess[0] = 4.0; a_guess[1] = 2.0; a_guess[2] = 3.0; 
	ia[3] = 1; ia[4] = 0; ia[5] = 0;
	a_guess[3] = 2.5; a_guess[4] = 5.0; a_guess[5] = 1.0;*/ 

	// Find one parameters
	/*ia[0] = 0; ia[1] = 1; ia[2] = 0;
	a_guess[0] = 5.0; a_guess[1] = 1.5; a_guess[2] = 3.0;
	ia[3] = 0; ia[4] = 0; ia[5] = 0;
	a_guess[3] = 3.0; a_guess[4] = 5.0; a_guess[5] = 1.0;*/

	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, fgauss, ITMAX, TOL, true); 

	// compute the residuals for the fit
	std::vector<std::vector<double>> data; 
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, fgauss, data); 

	// output the residuals 
	std::string thefile = "Gauss_non_lin_fit.txt"; 

	int nrows = 5; 
	vecut::write_into_file(thefile, data, nrows, npts);

	xdata.clear(); ydata.clear(); sigdata.clear(); 
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear(); 
	a.clear(); data.clear(); 
}

void testing::rng_test()
{
	// Examine the values of uniform random numbers generated by the RNG
	// Want to use this in the MC simulations
	// R. Sheehan 25 - 11 - 2019

	// Use current time information to generate seed value for rng
	time_t rawtime;

	time(&rawtime);

	int max, min, diff; 

	min = 7; max = 27; diff = max - min; // use these values to scale the RN to a certain range of integers

	long idum = rng::ranseed();
	// generate pairs of random numbers using the different rng
	for (int i = 0; i < 10; i++) {
		double RN = rng::ran1(&idum); 
		std::cout << min + diff * RN << " , " << std::round(min+diff*RN) << "\n";
	}
	std::cout << "\n"; 

	min = -7; max = -27; diff = max - min; // use these values to scale the RN to a certain range of integers

	idum = (-static_cast<long>(rawtime));
	// generate pairs of random numbers using the different rng
	for (int i = 0; i < 10; i++) {
		double RN = rng::ran1(&idum);
		std::cout << min + diff * RN << " , " << std::round(min + diff * RN) << "\n";
	}
	std::cout << "\n";

	max = 7; min = 0; 

	// generate pairs of random numbers using the different rng
	for (int i = 0; i < 11; i++) {
		std::cout << rng::ranint(&idum, min, max) << "\n";
	}
	std::cout << "\n";
}

void testing::monte_carlo_test()
{
	// Examine the monte-carlo method for estimating confidence limits on parameter estimates
	// R. Sheehan 25 - 11 - 2019

	// Use same example as before

	// Use current time information to generate seed value for rng
	time_t rawtime;

	time(&rawtime);

	// Generate data to use in the fit process
	int npts, npars;
	long idum = (-static_cast<long>(rawtime));
	double spread, xlow, xhigh, deltax, xpos, yval;

	xlow = 0.0; xhigh = 10.0; deltax = 0.1;
	npts = (int)(1.0 + ((xhigh - xlow) / deltax));

	std::vector<double> xdata(npts, 0.0);
	std::vector<double> ydata(npts, 0.0);
	std::vector<double> sigdata(npts, 0.0);

	npars = 6;
	std::vector<double> a(npars, 0.0);

	std::vector<double> dyda(npars, 0.0);

	// data stored in array a in the form a[i] = amp_{k}, a[i+1] = centre_{k}, a[i+2] = width_{k}
	a[0] = 5.0; a[1] = 2.0; a[2] = 3.0; // Gaussian 1
	a[3] = 3.0; a[4] = 5.0; a[5] = 1.0; // Gaussian 2

	spread = 0.1;
	xpos = xlow;
	for (int i = 0; i < npts; i++) {

		fgauss(xpos, a, &yval, dyda, npars); // evaluate the Gaussian function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		sigdata[i] = spread * yval;

		xpos += deltax;
	}

	a.clear(); // Don't need this anymore

	// Compute the estimate of the best fit to the data
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<double> a_guess_again(npars, 0.0); 
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	// Good guesses
	a_guess[0] = 4.5; a_guess[1] = 2.2; a_guess[2] = 2.8;
	a_guess[3] = 2.5; a_guess[4] = 4.3; a_guess[5] = 1.4;

	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, fgauss, ITMAX, TOL, true);

	/*xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();*/

	// Generate a random sample of the underlying data
	int nsmpl = static_cast<int>(0.9*npts); 
	std::vector<double> xs(nsmpl, 0.0);
	std::vector<double> ys(nsmpl, 0.0);
	std::vector<double> sigs(nsmpl, 0.0);
	std::vector<double> params(npars, 0.0);
	std::vector<double> errors(npars, 0.0);

	int ntest = 11; 

	for (int k = 0; k < ntest; k++) {

		fit::random_sample(xdata, ydata, sigdata, npts, xs, ys, sigs, nsmpl);

		// Perform a fit on the sampled data
		a_guess_again[0] = 4.5; a_guess_again[1] = 2.2; a_guess_again[2] = 2.8;
		a_guess_again[3] = 2.5; a_guess_again[4] = 4.3; a_guess_again[5] = 1.4;

		fit::non_lin_fit(xs, ys, sigs, nsmpl, a_guess_again, ia, npars, covar, alpha, &chisq, fgauss, ITMAX, TOL);

		for (int m = 0; m < npars; m++) {
			//std::cout << "da[" << m << "] = " << a_guess[m] - a_guess_again[m] << "\n";
			params[m] += a_guess_again[m];
			errors[m] += (a_guess[m] - a_guess_again[m]);
		}
		//std::cout << "\n";
	}

	std::cout << "Averaged Parameters with Errors from Monte-Carlo Sampling\n"; 
	for (int j = 0; j < npars; j++) {
		std::cout << "a[" << j << "] = " << params[j]/(ntest) << " +/- " << fabs(errors[j]/(ntest)) << "\n";
	}
}

void testing::diode_voltage(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// function that computes the diode voltage for input current x = I [mA]
	// a stores diode parameters a = { eta, T, Is}
	// a[0] = eta, a[1] = T, a[2] = Is
	// eta is the diode ideality, dimensionless constant between 1 - 10
	// T is the diode temperature, must be expressed in units of Celcius, converted to Kelvin during fit
	// Is diode saturation current, typically Is ~ 10^{-14} A = 10^{-11} mA
	// diode voltage value is given by *y
	// dyda is array that stores value of derivative of diode voltage function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 19 - 10 - 2021

	try {
		double kB_q = 8.61733e-5; // ( k_{B} / q ) in units of J / C K
		double Tkelvin = 273.15 + a[1]; // convert temperature to Kelvin and multiply by ( k_{B} / q ) 
		double Tterm = kB_q * Tkelvin; // convert temperature to Kelvin and multiply by ( k_{B} / q )
		//double Tterm = kB_q * a[1]; // convert temperature to Kelvin and multiply by ( k_{B} / q )
		double Iratio = (x / a[2]); // I_{d} / I_{s}
		double arg = 1.0 + Iratio; // compute term 1 + (I_{d} / I_{s})
		double log_term = log(arg); // log ( 1 + (I_{d} / I_{s}) )
		double mid_term = Tterm * log_term; // ( kB_q * Tkelvin ) * log ( 1 + (I_{d} / I_{s}) )
		*y = a[0] * mid_term; // eta * Tterm * log ( 1 + I_{d} / I_{s} )
		dyda[0] = mid_term; // \partial V_{d} / \partial \eta
		dyda[1] = a[0] * kB_q * log_term; // \partial V_{d} / \partial T
		dyda[2] = (-1.0 * a[0] * Tterm * Iratio) / (a[2] * arg); // \partial V_{d} / I_{s}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::resistive_diode_voltage(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// function that computes the diode voltage for input current x = I [mA] includes the contribution of a series resistance Rs [kOhm]
	// a stores diode parameters a = { eta, T, Is, Rs }
	// a[0] = eta, a[1] = T, a[2] = Is, a[3] = Rs
	// eta is the diode ideality, dimensionless constant between 1 - 10
	// T is the diode temperature, must be expressed in units of Celcius, converted to Kelvin during fit
	// Is diode saturation current, typically Is ~ 10^{-14} A = 10^{-11} mA
	// Rs is a resistor in series with the diode, expressed in units of kOhm, it's a convenient way to model the effect of metal contacts
	// diode voltage value is given by *y
	// dyda is array that stores value of derivative of diode voltage function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 7 - 4 - 2022

	try {
		double kB_q = 8.61733e-5; // ( k_{B} / q ) in units of J / C K
		double Tkelvin = 273.15 + a[1]; // convert temperature to Kelvin and multiply by ( k_{B} / q ) 
		double Tterm = kB_q * Tkelvin; // convert temperature to Kelvin and multiply by ( k_{B} / q )
		double Iratio = (x / a[2]); // I_{d} / I_{s}
		double arg = 1.0 + Iratio; // compute term 1 + (I_{d} / I_{s})
		double Rterm = x * a[3]; // R_{s} I_{d}
		double log_term = log(arg); // log ( 1 + (I_{d} / I_{s}) )
		double mid_term = Tterm * log_term; // ( kB_q * Tkelvin ) * log ( 1 + (I_{d} / I_{s}) )
		*y = a[0] * mid_term + Rterm; // eta * Tterm * log ( 1 + I_{d} / I_{s} ) + R_{s} I_{d}
		dyda[0] = mid_term; // \partial V_{d} / \partial \eta
		dyda[1] = a[0] * kB_q * log_term; // \partial V_{d} / \partial T
		dyda[2] = (-1.0 * a[0] * Tterm * Iratio) / (a[2] * arg); // \partial V_{d} / I_{s}
		dyda[3] = x; // \partial V_{d} / R_{s}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

//void testing::resistive_diode_voltage(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
//{
//	// function that computes the diode voltage for input current x = I [mA] includes the contribution of a series resistance Rs [kOhm]
//	// this is a simplified version, mostly interested in attempting to fit to the series resistance term
//	// a stores diode parameters a = { c1, c2, Rs }
//	// a[0] = c1 = (eta k_{B} T / q), a[1] = c2 = I_{s}, a[2] = Rs
//	// eta is the diode ideality, dimensionless constant between 1 - 10
//	// T is the diode temperature, must be expressed in units of Celcius, converted to Kelvin during fit
//	// Is diode saturation current, typically Is ~ 10^{-14} A = 10^{-11} mA
//	// Rs is a resistor in series with the diode, expressed in units of kOhm, it's a convenient way to model the effect of metal contacts
//	// diode voltage value is given by *y
//	// dyda is array that stores value of derivative of diode voltage function wrt each parameter in a
//	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
//	// na is no. parameters
//	// R. Sheehan 7 - 4 - 2022
//
//	try {
//		double Iratio = (x / a[1]); // I_{d} / c2
//		double arg = 1.0 + Iratio; // 1 + ( I_{d} / c2 )
//		double log_arg = log(arg); // log ( 1 + ( I_{d} / c2 ) )
//		*y = ( a[0] * log_arg ) + ( a[2] * x ); // c1 * log ( 1 + I_{d} / c2 ) + R_{s} I_{d}
//		dyda[0] = log_arg; // \partial V_{d} / \partial c1
//		dyda[1] = ( -1.0 * ( a[0] / a[1] ) * Iratio )/arg; // \partial V_{d} / \partial c2		
//		dyda[2] = x; // \partial V_{d} / R_{s}
//	}
//	catch (std::invalid_argument& e) {
//		std::cerr << e.what();
//	}
//}

void testing::diode_voltage_data_fit()
{
	// Apply LM method to measured diode voltage data
	// R. Sheehan 19 - 10 - 2021

	double eta = 1.5; // diode ideality (unitless)
	double T = 22; // diode temperature ( Kelvin )
	double Is = 5e-11; // diode saturation current (mA)
	double Id = 50; // diode current (mA)
	double y = 0.0; //computed voltage value

	int npars = 3;
	std::vector<double> a(npars, 0.0);
	std::vector<double> dyda(npars, 0.0);

	// Initial guesses for the parameters
	a[0] = eta; a[1] = T; a[2] = Is;

	testing::diode_voltage(Id, a, &y, dyda, npars);

	std::cout << "Diode Fit Test Evaluation\n"; 
	std::cout << "Diode current value: " << Id << " mA\n"; 
	std::cout << "Diode voltage value: " << y << " V\n"; 
	std::cout << "Der wrt eta: " << dyda[0] << "\n"; 
	std::cout << "Der wrt T: " << dyda[1] << "\n"; 
	std::cout << "Der wrt Is: " << dyda[2] << "\n\n"; 
	std::cout << "Diode Equation Fit\n"; 

	// Generate data to use in the fit process
	int npts;
	long idum = (-1011);
	double spread, xlow, xhigh, deltax, xpos, yval;

	// Diode fit works better when you ignore the value at 0 mA
	// 
	// You've noticed this before
	xlow = 0.45; xhigh = 100.0; deltax = 0.45;
	npts = (int)(1.0 + ((xhigh - xlow) / deltax));

	std::vector<double> xdata(npts, 0.0);
	std::vector<double> ydata(npts, 0.0);
	std::vector<double> sigdata(npts, 0.0);

	spread = 0.001; // variance of the noise being added to the signal
	xpos = xlow;
	for (int i = 0; i < npts; i++) {

		diode_voltage(xpos, a, &yval, dyda, npars); // evaluate the diode-voltage function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		sigdata[i] = fabs(yval) > 0.0 ? spread * yval : spread; // sigdata cannot have zero values

		xpos += deltax;
	}

	// Perform the best it search for the data set
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	ia[1] = 0; // search for params 0 and 2, fix param 1 value
	a_guess[0] = 10; a_guess[1] = 22.0; a_guess[2] = 4e-10; // initial guesses for the parameters

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, diode_voltage, ITMAX, TOL, true);

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, diode_voltage, data);

	// output the residuals 
	std::string thefile = "Diode_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	a.clear();
}

void testing::diode_voltage_data_fit_test()
{
	// Apply the LM method to measured diode voltage data
	// R. Sheehan 15 - 12 - 2021

	// Read in the measured diode voltage data
	//std::string filename = "AM_FR101_Swp_3_Rd_10.txt";
	std::string filename = "AM_IN4001_Swp_3_Rd_10.txt";

	int npts, n_rows, npars = 3, n_cols, indx_max = 0;
	double spread = 0.05, I_start = 0.5;

	time_t rawtime;
	time(&rawtime);
	long idum = (-static_cast<long>(rawtime)); // randomised seed for rng

	std::vector<std::vector<double>> the_data;

	vecut::read_into_matrix(filename, the_data, n_rows, n_cols, true);

	// Estimate sigdata
	std::vector<double> xdata;
	std::vector<double> ydata;
	std::vector<double> sigdata;

	// for some reason the copy was throwing an exception
	//std::copy( xdata.begin(), xdata.end(), vecut::get_col(the_data, 0) ); 
	//std::copy( ydata.begin(), ydata.end(), vecut::get_col(the_data, 1) ); 

	//xdata = vecut::get_col(the_data, 0);
	//ydata = vecut::get_col(the_data, 1);

	for (int i = 0; i < n_rows; i++) {
		if (the_data[i][0] > I_start) {
			xdata.push_back(the_data[i][0]); // read in current data for I > 0 mA
			ydata.push_back(the_data[i][1]); // store voltage data
			if (fabs(the_data[i][1]) > 0.0) {
				sigdata.push_back(spread * the_data[i][1]); // estimate error in voltage reading
			}
			else {
				sigdata.push_back(spread); // sigdata cannot have zero values
			}
		}
	}

	npts = static_cast<int>(xdata.size());

	// Perform the best it search for the data set
	int ITMAX = 50;

	double TOL = 1e-3;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	ia[1] = 0; 
	a_guess[0] = 10.0; a_guess[1] = 20.0; a_guess[2] = 1e-5; // initial guesses for the parameters

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, diode_voltage, ITMAX, TOL, true);

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, diode_voltage, data);

	// output the residuals 
	std::string thefile = "Diode_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	// take a look at the goodness of fit statistics
	double chisqr = 0.0, rsqr = 0.0, dof = static_cast<int>(npts - npars), gof = 0.0;
	fit::goodness_of_fit(xdata, ydata, data[4], npts, a_guess, npars, diode_voltage, &chisqr, &dof, &rsqr, &gof);

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	data.clear(); the_data.clear();
}

void testing::resistive_diode_voltage_data_fit()
{
	// Apply LM method to measured diode voltage data
	// R. Sheehan 7 - 4 - 2022

	double eta = 1.5; // diode ideality (unitless)
	double T = 22; // diode temperature ( Celcius )
	double Is = 5e-11; // diode saturation current (mA)
	double Id = 50; // diode current (mA)
	double Rs = (3.7 / 1000.0); // diode series resistance (kOhm)
	double y = 0.0; //computed voltage value

	int npars = 4;
	std::vector<double> a(npars, 0.0);
	std::vector<double> dyda(npars, 0.0);

	// Initial guesses for the parameters
	a[0] = eta; a[1] = T; a[2] = Is; a[3] = Rs; 

	testing::resistive_diode_voltage(Id, a, &y, dyda, npars);

	std::cout << "Diode Fit Test Evaluation\n";
	std::cout << "Diode current value: " << Id << " mA\n";
	std::cout << "Diode voltage value: " << y << " V\n";
	std::cout << "Der wrt eta: " << dyda[0] << "\n";
	std::cout << "Der wrt T: " << dyda[1] << "\n";
	std::cout << "Der wrt Is: " << dyda[2] << "\n";
	std::cout << "Der wrt Rs: " << dyda[3] << "\n\n";
	std::cout << "Diode Equation Fit\n";

	// Generate data to use in the fit process
	int npts;
	long idum = (-1011);
	double spread, xlow, xhigh, deltax, xpos, yval;

	// Diode fit works better when you ignore the value at 0 mA
	// 
	// You've noticed this before
	xlow = 0.5; xhigh = 100.0; deltax = 0.45;
	npts = (int)(1.0 + ((xhigh - xlow) / deltax));

	std::vector<double> xdata(npts, 0.0);
	std::vector<double> ydata(npts, 0.0);
	std::vector<double> sigdata(npts, 0.0);

	spread = 0.01; // variance of the noise being added to the signal
	xpos = xlow;
	for (int i = 0; i < npts; i++) {

		resistive_diode_voltage(xpos, a, &yval, dyda, npars); // evaluate the diode-voltage function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		sigdata[i] = fabs(yval) > 0.0 ? spread * yval : spread; // sigdata cannot have zero values

		xpos += deltax;
	}

	// Perform the best it search for the data set
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	ia[1] = 0; // search for params 0 and 2, fix param 1 value
	a_guess[0] = 10; a_guess[1] = T; a_guess[2] = 4e-10; a_guess[3] = (20.0/1000.0); // initial guesses for the parameters

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, resistive_diode_voltage, ITMAX, TOL, true);

	std::cout << "\nEstimate of the series resistance: " << 1000.0 * a_guess[3] << " Ohm\n\n"; 

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, resistive_diode_voltage, data);

	// output the residuals 
	std::string thefile = "Diode_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	a.clear();
}

void testing::Lorentzian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// Definition of the Lorentzian function to be fitted
	// a stores Lorentzian parameters a = { A, x_{centre}, G/2}
	// a[0] = A, a[1] = x_{centre}, a[2] = G / 2
	// G/2 is the Lorentzian half-width at half-maximum (i.e. linewidth)
	// Lorentzian value is given by *y
	// This Lorentzian is not normalised, to normalise multiply *y by 1/pi
	// Normalisation not required for my purposes
	// dyda is array that stores value of derivative of Lorentzian function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 21 - 10 - 2021

	try {
		double t1 = x - a[1]; // ( x - x_{centre} )
		double t2 = (a[0] * a[2]); // A (G/2)
		double denom = ( template_funcs::DSQR(t1) + template_funcs::DSQR(a[2]) ); // ( x - x_{centre} )^{2} + (G/2)^{2}
		*y = t2 / denom; // A (G/2) / [ ( x - x_{centre} )^{2} + (G/2)^{2} ]
		dyda[0] = (*y) / a[0]; // \partial L / \partial A
		dyda[1] = (2.0 * t1 * template_funcs::DSQR(*y) ) / t1; // \partial L / \partial x_{centre}
		dyda[2] = (*y) * ( ( 1.0 / a[2]) -  ( ( 2.0 * (*y) ) / a[0] ) ); // \partial L / \partial G
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::Lorentzian_data_fit()
{
	// Apply LM method to test Lorentzian data
	// R. Sheehan 21 - 10 - 2021

	double f = 80; 
	double A = 1; // Lorentzian amplitude
	double xc = 80; // Lorentzian centre
	double G2 = 0.74; // Lorentzian HWHM
	double y = 0.0; //computed Lorentzian value

	int npars = 3;
	std::vector<double> a(npars, 0.0);
	std::vector<double> dyda(npars, 0.0);

	// Initial guesses for the parameters
	a[0] = A; a[1] = xc; a[2] = G2;

	/*testing::Lorentzian(f, a, &y, dyda, npars);

	std::cout << "Lorentzian Fit Test Evaluation\n";
	std::cout << "Lorentzian Frequency value: " << f << " MHz\n";
	std::cout << "Lorentizian Spectrum value: " << y << "\n";
	std::cout << "Der wrt xc: " << dyda[0] << "\n";
	std::cout << "Der wrt G: " << dyda[1] << "\n\n";*/
	std::cout << "Lorentzian Equation Fit\n";

	// Generate data to use in the fit process
	int npts;
	double spread, xlow, xhigh, deltax, xpos, yval;

	time_t rawtime;
	time(&rawtime);
	long idum = (-static_cast<long>(rawtime)); // randomised seed for rng

	xlow = 70.0; xhigh = 90.0; deltax = 0.1;
	npts = (int)(1.0 + ((xhigh - xlow) / deltax));

	std::vector<double> xdata(npts, 0.0);
	std::vector<double> ydata(npts, 0.0);
	std::vector<double> sigdata(npts, 0.0);

	spread = 0.1; // variance of the noise being added to the signal
	xpos = xlow;
	for (int i = 0; i < npts; i++) {

		Lorentzian(xpos, a, &yval, dyda, npars); // evaluate the Lorentzian function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		sigdata[i] = fabs(yval) > 0.0 ? spread * yval : spread; // sigdata cannot have zero values
		//sigdata[i] = spread/100.0; // sigdata cannot have zero values
		//sigdata[i] = 1.0; // sigdata cannot have zero values

		xpos += deltax;
	}

	// For the Lorentzian the HWHM value is given by 1/peak-value
	// So for noisy data HWHM should be roughly equal to 1/peak-value of data
	// Max value should be close to x_centre, assuming data is distributed equally around x_centre
	// This doesn't work when x_centre is outside the range [xlow, xhigh], but otherwise works quite well
	int imax=0; 
	double Lmax = -500; 
	for (int i = 0; i < npts; i++) {
		if (ydata[i] > Lmax) {
			Lmax = ydata[i];
			imax = i; 
		}
	}
	std::cout << "Max value in data set: " << Lmax << "\n"; 
	std::cout << "Corresponding Frequency: " << xdata[imax] << "\n"; 
	std::cout << "Estimate of HWHM: " << 1.0 / Lmax << "\n\n"; 

	// Perform the best it search for the data set
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	//ia[0] = 0; // search for params 0 and 2, fix param 1 value
	//a_guess[0] = A; a_guess[1] = xc; a_guess[2] = 4.0; // initial guesses for the parameters, fit routine is sufficiently robust
	a_guess[0] = A; a_guess[1] = xc; a_guess[2] = 1.0 / Lmax; // initial guesses for the parameters

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, Lorentzian, ITMAX, TOL, true);

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, Lorentzian, data);

	// output the residuals 
	std::string thefile = "Lorentzian_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	// take a look at the goodness of fit statistics
	double chisqr = 0.0, rsqr = 0.0, dof = static_cast<int>(npts - npars), gof = 0.0; 
	fit::goodness_of_fit(xdata, ydata, sigdata, npts, a_guess, npars, Lorentzian, &chisqr, &dof, &rsqr, &gof); 

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	a.clear(); data.clear(); 
}

void testing::Lorentzian_data_fit_test()
{
	// Apply LM method to measured linewidth spectrum data
	// R. Sheehan 26 - 10 - 2021

	// Read in the measured spectral data
	//std::string filename = "Sample_LLM.csv"; 
	//std::string filename = "Lorentz_iodeal.csv"; // this is the same data set as Sample_LLM.csv
	//std::string filename = "Smpl_LLM_7.txt";
	//std::string filename = "LLM_Spctrm_I_65.txt";

	/*std::string dir_name = "C:\\users\\robertsheehan\\Research\\Laser_Physics\\Linewidth\\Data\\ESA_Spectra_Versus_VOA_Bias\\";
	useful_funcs::set_directory(dir_name);
	useful_funcs::get_directory();

	std::string filename = "JDSU_DFB_T_20_I_50_V_300.txt";*/

	std::string dir_name = "C:\\users\\robertsheehan\\Research\\Laser_Physics\\Linewidth\\Data\\NKT_LCR_DSHI\\";
	useful_funcs::set_directory(dir_name);
	useful_funcs::get_directory();

	bool ReScaleFrq = true;
	std::string unitsFrq = ReScaleFrq ? " kHz" : " MHz";
	std::string unitsPow = " pW";

	int fbeat = 320; // beat frequency MHz
	std::string filename = "NKT_I_100_Vb_30_RBW_05_fb_" + template_funcs::toString(fbeat) + dottxt;

	int npts, n_rows, npars = 3, n_cols, indx_max = 0;
	long idum = (-1011);
	double spread = 0.15, spctr_max = -500.0, f_max = 0, f_start, f_end, scale_fac;
	
	std::vector<std::vector<double>> the_data; 

	vecut::read_into_matrix(filename, the_data, n_rows, n_cols, true);

	// Estimate sigdata
	std::vector<double> xdata;
	std::vector<double> ydata;
	std::vector<double> sigdata;

	// for some reason the copy was throwing an exception
	//std::copy( xdata.begin(), xdata.end(), vecut::get_col(the_data, 0) ); 
	//std::copy( ydata.begin(), ydata.end(), vecut::get_col(the_data, 1) ); 

	//xdata = vecut::get_col(the_data, 0);
	//ydata = vecut::get_col(the_data, 1);

	scale_fac = 1.0e+6; // why is this scale factor being deployed? To convert power from dBm to pW
	//scale_fac = 1.0e+3; // why is this scale factor being deployed? To convert power from dBm to nW
	for (int i = 0; i < n_rows; i++) {
		if (ReScaleFrq) {
			xdata.push_back(1000 * (the_data[i][0] - (double)(fbeat))); // rescale frq to zero in units of kHz
		}
		else {
			xdata.push_back(the_data[i][0]); // rescale frq to zero in units of kHz			
		}
		ydata.push_back(scale_fac * convert_dBm_to_mW(the_data[i][1])); // convert the spectral data from dBm to mW scale
	}

	npts = static_cast<int>( xdata.size() ); 

	for (int i = 0; i < npts; i++) {
		if (ydata[i] > spctr_max) {
			spctr_max = ydata[i]; 
			indx_max = i; 
		}
		
		if (fabs(ydata[i]) > 0.0) {
			sigdata.push_back(spread * ydata[i]); 
		}
		else {
			sigdata.push_back(spread); // sigdata cannot have zero values
		}
	}

	// Scale peak to unity
	// This won't work because it will force the HWHM to be equal to 0.5 or something
	/*for (int i = 0; i < npts; i++) {
		ydata[i] /= spctr_max; 
	}*/

	// What about subtracting an offset? Doesn't work either
	/*for (int i = 0; i < npts; i++) {
		ydata[i] -= spctr_max; 
	}*/

	// For the Lorentzian the HWHM value is given by 1/peak-value
	// So for noisy data HWHM should be roughly equal to 1/peak-value of data
	// Max value should be close to x_centre, assuming data is distributed equally around x_centre
	// This doesn't work when x_centre is outside the range [xlow, xhigh], but otherwise works quite well
	
	f_max = xdata[indx_max]; // find location of peak value

	std::cout << "\nMax value in data set: " << spctr_max << unitsPow << "\n";
	std::cout << "Corresponding Frequency: " << f_max << unitsFrq << "\n\n";

	// Perform the best it search for the data set
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	// initial guesses for the parameters	switch to tell code whether or not they're to be fitted
	//a_guess[0] = ReScaleFrq ? spctr_max : 10 * spctr_max;			ia[0] = ReScaleFrq ? 0 : 1;
	a_guess[0] = 10 * spctr_max;			ia[0] = 1;
	a_guess[1] = f_max;				ia[1] = 0;
	a_guess[2] = 3.32;				ia[2] = 1;

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, Lorentzian, ITMAX, TOL, true);

	std::cout << "Fitted centre freq: " << a_guess[1] << unitsFrq << "\n";
	std::cout << "Computed peak val: " << a_guess[0] / a_guess[2] << unitsPow << "\n";
	std::cout << "Lorentz HWHM: " << a_guess[2] << unitsFrq << "\n";

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, Lorentzian, data);

	// output the residuals 
	std::string thefile = "Lorentzian_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	// take a look at the goodness of fit statistics
	double chisqr = 0.0, rsqr = 0.0, dof = static_cast<int>(npts - npars), gof = 0.0;
	fit::goodness_of_fit(xdata, ydata, data[4], npts, a_guess, npars, Lorentzian, &chisqr, &dof, &rsqr, &gof);
}

void testing::Gaussian(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// Definition of the Gaussian function to be fitted
	// a stores Gaussian parameters a = { amplitude, mean, standard deviation} = {A, b, c}
	// a[0] = A, a[1] = b, a[2] = c
	// c is related the Gaussian half-width at half-maximum (i.e. linewidth) HWHM = sqrt{ 2 log(2) } c ~ 1.17741 c
	// Gaussian value is given by *y
	// This Gaussian is not normalised, to normalise multiply *y by 1/(a c sqrt{pi})
	// Normalisation not required for my purposes
	// dyda is array that stores value of derivative of Lorentzian function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 27 - 10 - 2021

	try {
		double t1 = x - a[1]; // ( x - b )
		double t1sqr = template_funcs::DSQR(t1); // ( x - b )^{2}
		double csqr = template_funcs::DSQR(a[2]); // c^{2}
		double arg = (-1.0 * t1sqr) / (2.0*csqr); // -( x - b )^{2} / 2 c^{2}
		*y = a[0] * exp(arg); // A exp( -( x - b )^{2} / 2 c^{2} )
		dyda[0] = (*y)/a[0]; // \partial G / \partial A
		dyda[1] = (t1 / csqr) * (*y); // \partial G / \partial b
		dyda[2] = (t1sqr/(a[2]*csqr))*(*y); // \partial G / \partial c
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::Gaussian_data_fit()
{
	// Apply LM method to test Gaussian data
	// R. Sheehan 27 - 10 - 2021

	double f = 82;
	double A = 1.5; // Amplitude
	double b = 80; // Gaussian centre
	double c = 0.5; // Gaussian width
	double y = 0.0; //computed Gaussian value

	int npars = 3;
	std::vector<double> a(npars, 0.0);
	std::vector<double> dyda(npars, 0.0);

	// Initial guesses for the parameters
	a[0] = A; a[1] = b; a[2] = c; 

	testing::Gaussian(f, a, &y, dyda, npars);

	std::cout << "Gaussian Fit Test Evaluation\n";
	std::cout << "Gaussian Frequency value: " << f << " MHz\n";
	std::cout << "Gaussian Spectrum value: " << y << "\n";
	std::cout << "Der wrt A: " << dyda[0] << "\n";
	std::cout << "Der wrt b: " << dyda[1] << "\n";
	std::cout << "Der wrt c: " << dyda[2] << "\n\n";
	std::cout << "Gaussian Equation Fit\n";

	// Generate data to use in the fit process
	int npts;
	double spread, xlow, xhigh, deltax, xpos, yval;

	time_t rawtime;
	time(&rawtime);
	long idum = (-static_cast<long>(rawtime)); // randomised seed for rng

	xlow = 75.0; xhigh = 85.0; deltax = 0.1;
	npts = (int)(1.0 + ((xhigh - xlow) / deltax));

	std::vector<double> xdata(npts, 0.0);
	std::vector<double> ydata(npts, 0.0);
	std::vector<double> sigdata(npts, 0.0);

	spread = 0.001; // variance of the noise being added to the signal
	xpos = xlow;
	for (int i = 0; i < npts; i++) {

		Gaussian(xpos, a, &yval, dyda, npars); // evaluate the Lorentzian function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		sigdata[i] = fabs(yval) > 0.0 ? spread * yval : spread; // sigdata cannot have zero values

		xpos += deltax;
	}

	// Perform the best it search for the data set
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	ia[0] = 0; ia[1] = 0; // search for params 0 and 2, fix param 1 value
	a_guess[0] = A; a_guess[1] = b; a_guess[2] = 1.0; // initial guesses for the parameters, fit routine is sufficiently robust

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, Gaussian, ITMAX, TOL, true);

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, Gaussian, data);

	// output the residuals 
	std::string thefile = "Gaussian_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	a.clear(); data.clear();
}

void testing::Gaussian_data_fit_test()
{
	// Apply LM method to measured linewidth spectrum data
	// R. Sheehan 27 - 10 - 2021

	// Read in the measured spectral data
	std::string filename = "Smpl_LLM_1.txt";
	//std::string filename = "Smpl_LLM_2.txt";

	int npts, n_rows, npars = 3, n_cols, indx_max = 0;
	long idum = (-1011);
	double spread = 0.01, spctr_max = -500.0, f_max = 0, f_start, f_end, scale_fac;

	std::vector<std::vector<double>> the_data;

	vecut::read_into_matrix(filename, the_data, n_rows, n_cols, true);

	// Estimate sigdata
	std::vector<double> xdata;
	std::vector<double> ydata;
	std::vector<double> sigdata;

	// for some reason the copy was throwing an exception
	//std::copy( xdata.begin(), xdata.end(), vecut::get_col(the_data, 0) ); 
	//std::copy( ydata.begin(), ydata.end(), vecut::get_col(the_data, 1) ); 

	//xdata = vecut::get_col(the_data, 0);
	//ydata = vecut::get_col(the_data, 1);

	scale_fac = 1.0e+6;  f_start = 70.0; f_end = 90.0;
	for (int i = 0; i < n_rows; i++) {
		if (the_data[i][0] > f_start && the_data[i][0] < f_end) {
			xdata.push_back(the_data[i][0]);
			ydata.push_back(scale_fac * pow(10.0, the_data[i][1] / 10.0)); // convert the spectral data from dBm to mW scale and rescale it
		}
	}

	npts = static_cast<int>(xdata.size());

	for (int i = 0; i < npts; i++) {
		if (ydata[i] > spctr_max) {
			spctr_max = ydata[i];
			indx_max = i;
		}

		if (fabs(ydata[i]) > 0.0) {
			sigdata.push_back(spread * ydata[i]);
		}
		else {
			sigdata.push_back(spread); // sigdata cannot have zero values
		}
	}

	std::cout << "Max value in data set: " << spctr_max << "\n";
	std::cout << "Corresponding Frequency: " << xdata[indx_max] << "\n\n";

	// Perform the best it search for the data set
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	ia[0] = 0; ia[1] = 0;
	a_guess[0] = spctr_max;
	//a_guess[1] = xdata[indx_max];
	a_guess[1] = 80;
	a_guess[2] = 1.0; // initial guesses for the parameters

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, Gaussian, ITMAX, TOL, true);

	std::cout << "FWHM: " << 2.0 * sqrt(2.0 * log(2.0)) * a_guess[2] << "\n"; 
	std::cout << "HWHM: " << sqrt(2.0 * log(2.0)) * a_guess[2] << "\n\n"; 

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, Gaussian, data);

	// output the residuals 
	std::string thefile = "Gaussian_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	// take a look at the goodness of fit statistics
	double chisqr = 0.0, rsqr = 0.0, dof = static_cast<int>(npts - npars), gof = 0.0;
	fit::goodness_of_fit(xdata, ydata, data[4], npts, a_guess, npars, Gaussian, &chisqr, &dof, &rsqr, &gof);

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	data.clear(); the_data.clear();
}

void testing::Voigt(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// Definition of the Voigt function to be fitted
	// a stores Voigt parameters a = { h, x_{centre}, g, sigma}
	// a[0] = h, a[1] = x_{centre}, a[2] = g, a[3] = sigma
	// h is an amplitude fitting factor
	// x_centre is the centre of the Voigt frunction
	// g is the half-width at half-maximum of the Lorentzian portion of Voigt
	// sigma is the std. dev. of the Gaussian portion of Voigt HWHM_{Gauss} = sqrt( 2 log(2) ) c
	// It should be possible to express HWHM_{Voigt} in terms of g and sigma
	// Voigt value is given by *y
	// dyda is array that stores value of derivative of Voigt function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 30 - 11 - 2021

	try {
		std::complex<double> z = ( (x - a[1] + eye * a[2]) / a[3]  ); // z = (x - x_{0} + i g) / sigma
		std::complex<double> W = Faddeeva::w(z); // w(z) = exp(-z^{2}) Erfc(-i z)
		std::complex<double> dW = ((2.0 * eye) / SQRT_PI) - ( 2.0 * z * W ); // dw / dz by definition
		double h_sig = a[0] / a[3]; // h / sigma

		*y = real( a[0] * W ); // V = re( h w(z) )
		dyda[0] = real(W); // \partial V / \partial h
		dyda[1] = real(-1.0 * h_sig * dW); // \partial V / \partial x_{0} = - \partial V / \partial x
		dyda[2] = real(eye * h_sig * dW); // \partial V / \partial g
		dyda[3] = real(-1.0 * h_sig * z * dW); // \partial V / \partial sigma
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::Voigt_data_fit()
{
	// Apply LM method to test Voigt data
	// R. Sheehan 30 - 11 - 2021

	double f = 81;
	double h = 1.5; // Amplitude
	double x0 = 80; // Voigt centre
	double g = 1.0; // Lorentzian HWHM
	double sigma = 1.5; // Gaussian std. dev.
	double y = 0.0; //computed Voigt value

	int npars = 4;
	std::vector<double> a(npars, 0.0);
	std::vector<double> dyda(npars, 0.0);

	// Initial guesses for the parameters
	a[0] = h; a[1] = x0; a[2] = g; a[3] = sigma;

	testing::Voigt(f, a, &y, dyda, npars);

	std::cout << "Voigt Fit Test Evaluation\n";
	std::cout << "Voigt Frequency value: " << f << " MHz\n";
	std::cout << "Voigt Spectrum value: " << y << "\n";
	std::cout << "Der wrt h: " << dyda[0] << "\n";
	std::cout << "Der wrt x0: " << dyda[1] << "\n";
	std::cout << "Der wrt g: " << dyda[2] << "\n";
	std::cout << "Der wrt sigma: " << dyda[3] << "\n\n";
	std::cout << "Voigt Equation Fit\n";

	// Generate data to use in the fit process
	int npts;
	double spread, xlow, xhigh, deltax, xpos, yval;

	time_t rawtime;
	time(&rawtime);
	long idum = (-static_cast<long>(rawtime)); // randomised seed for rng

	xlow = 70.0; xhigh = 90.0; deltax = 0.01;
	npts = (int)(1.0 + ((xhigh - xlow) / deltax));

	std::vector<double> xdata(npts, 0.0);
	std::vector<double> ydata(npts, 0.0);
	std::vector<double> sigdata(npts, 0.0);

	spread = 0.1; // variance of the noise being added to the signal
	xpos = xlow;
	for (int i = 0; i < npts; i++) {

		Voigt(xpos, a, &yval, dyda, npars); // evaluate the Lorentzian function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		sigdata[i] = fabs(yval) > 0.0 ? spread * yval : spread; // sigdata cannot have zero values

		xpos += deltax;
	}

	// Perform the best it search for the data set
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	//ia[1] = 0; // search for params 0 and 2, fix param 1 value
	a_guess[0] = 3*h; 
	a_guess[1] = x0; 
	a_guess[2] = 2*g; 
	a_guess[3] = 2*g; // initial guesses for the parameters, fit routine is sufficiently robust

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, Voigt, ITMAX, TOL, true);

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, Voigt, data);

	// output the residuals 
	std::string thefile = "Voigt_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	a.clear(); data.clear();
}

void testing::Voigt_data_fit_test()
{
	// Apply LM method to measured linewidth spectrum data
	// 
	// R. Sheehan 27 - 10 - 2021

	// Read in the measured spectral data
	//std::string filename = "Smpl_LLM_11.txt";
	//std::string filename = "LLM_Spctrm_I_60.txt";

	//std::string dir_name = "C:\\users\\robertsheehan\\Research\\Laser_Physics\\Linewidth\\Data\\ESA_Spectra_Versus_VOA_Bias\\"; 
	std::string dir_name = "C:\\users\\robertsheehan\\Research\\Laser_Physics\\Linewidth\\Data\\NKT_LCR_DSHI\\"; 
	useful_funcs::set_directory(dir_name);
	useful_funcs::get_directory(); 

	//std::string filename = "JDSU_DFB_T_20_I_50_V_300.txt"; 

	bool ReScaleFrq = true; 
	std::string unitsFrq = ReScaleFrq ? " kHz" : " MHz";
	std::string unitsPow = " pW";

	int fbeat = 400; // beat frequency MHz
	std::string filename = "NKT_I_100_Vb_30_RBW_05_fb_" + template_funcs::toString(fbeat) + dottxt; 
	
	int npts, n_rows, npars = 4, n_cols, indx_max = 0;
	double spread = 0.05, spctr_max = -500.0, f_max = 0, scale_fac;

	time_t rawtime;
	time(&rawtime);
	long idum = (-static_cast<long>(rawtime)); // randomised seed for rng

	std::vector<std::vector<double>> the_data;

	vecut::read_into_matrix(filename, the_data, n_rows, n_cols, true);

	// Estimate sigdata
	std::vector<double> xdata;
	std::vector<double> ydata;
	std::vector<double> sigdata;

	// for some reason the copy was throwing an exception
	//std::copy( xdata.begin(), xdata.end(), vecut::get_col(the_data, 0) ); 
	//std::copy( ydata.begin(), ydata.end(), vecut::get_col(the_data, 1) ); 

	//xdata = vecut::get_col(the_data, 0);
	//ydata = vecut::get_col(the_data, 1);

	scale_fac = 1.0e+6; // why is this scale factor being deployed? To convert power from dBm to pW
	//scale_fac = 1.0e+3; // why is this scale factor being deployed? To convert power from dBm to nW
	for (int i = 0; i < n_rows; i++) {
		if (ReScaleFrq) {
			xdata.push_back( 1000 * (the_data[i][0] - (double)(fbeat)) ); // rescale frq to zero in units of kHz
		}
		else {
			xdata.push_back( the_data[i][0] ); // rescale frq to zero in units of kHz			
		}
		ydata.push_back(scale_fac * convert_dBm_to_mW(the_data[i][1])); // convert the spectral data from dBm to mW scale
	}

	npts = static_cast<int>(xdata.size());

	// find location of peak and its value
	for (int i = 0; i < npts; i++) {
		if (ydata[i] > spctr_max) {
			spctr_max = ydata[i];
			indx_max = i;
		}

		if (fabs(ydata[i]) > 0.0) {
			sigdata.push_back(spread * ydata[i]);
		}
		else {
			sigdata.push_back(spread); // sigdata cannot have zero values
		}
	}

	f_max = xdata[indx_max]; // find location of peak value

	std::cout << "\nMax value in data set: " << spctr_max << unitsPow << "\n";
	std::cout << "Corresponding Frequency: " << f_max << unitsFrq << "\n\n";

	// Perform the best it search for the data set
	int ITMAX = 50;

	double TOL = 1e-3;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	// initial guesses for the parameters	switch to tell code whether or not they're to be fitted
	a_guess[0] = ReScaleFrq ? spctr_max : 10*spctr_max;			ia[0] = ReScaleFrq ? 0 : 1;
	a_guess[1] = f_max;				ia[1] = 0;
	a_guess[2] = 3.32;				ia[2] = 1; 
	a_guess[3] = 3.32;				ia[3] = 1; 

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, Voigt, ITMAX, TOL, true);

	// compute the HWHM using bisection
	double xlow = a_guess[1], xhigh = xdata[npts-1], HWHM = 0.0;
	
	Voigt_HWHM(xlow, xhigh, a_guess, npars, &HWHM);

	std::cout << "Fitted centre freq: " << a_guess[1] << unitsFrq <<"\n";
	std::cout << "Computed peak val: " << a_guess[0] * exp( template_funcs::DSQR(a_guess[2]/ a_guess[3]) ) * probability::erffc(a_guess[2] / a_guess[3]) << unitsPow <<"\n";
	std::cout << "Lorentz HWHM: " << a_guess[2] << unitsFrq << "\n";
	std::cout << "Gauss HWHM: " << sqrt( 2.0*log(2.0) ) * a_guess[3] << unitsFrq <<"\n";
	std::cout << "Voigt HWHM: " << HWHM << unitsFrq << "\n\n";

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, Voigt, data);

	// output the residuals 
	std::string thefile = "Voigt_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	// take a look at the goodness of fit statistics
	double chisqr = 0.0, rsqr = 0.0, dof = static_cast<int>(npts - npars), gof = 0.0;
	fit::goodness_of_fit(xdata, ydata, data[4], npts, a_guess, npars, Voigt, &chisqr, &dof, &rsqr, &gof);

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	data.clear(); the_data.clear();
}

void testing::Voigt_HWHM(double xlow, double xhigh, std::vector<double>& a, int& na, double* HWHM, bool loud)
{
	// determine the HWHM of a Voigt function with the parameter set defined in a
	// use bisection method algorithm to look for the HWHM on the interval [xlow, xhigh]
	// a stores Voigt parameters a = { h, x_{centre}, g, sigma}
	// a[0] = h, a[1] = x_{centre}, a[2] = g, a[3] = sigma
	// h is an amplitude fitting factor
	// x_centre is the centre of the Voigt frunction
	// g is the half-width at half-maximum of the Lorentzian portion of Voigt
	// sigma is the std. dev. of the Gaussian portion of Voigt HWHM_{Gauss} = sqrt( 2 log(2) ) c
	// parameters in a will have been determined by Levenberg-Marquardt non-linear fit procedure
	// result will be stored in HWHM
	// R. Sheehan 1 - 12 - 2021

	try {
	
		bool c1 = xhigh > xlow ? true : false; 
		bool c2 = xlow > 0 ? true : false; 
		bool c3 = a[3] > 0 ? true : false; 
		bool c10 = c1 && c3;

		if (c10) {
			int count = 0, MAXIT = 50; // max. no iterations of bisection method
			bool converged = false; 
			double TOL = 1.0e-3; // desired accuracy of root computation
			double root = xlow, left = xlow, right = xhigh, dx = 0.0, fl = 0.0, fr = 0.0; 
			double lor_gau = a[2] / a[3]; // g_{lor} / sigma_{gau}
			double Vmax_half = 0.5*( a[0] * exp(template_funcs::DSQR(lor_gau)) * probability::erffc(lor_gau) ); // half the peak value of the Voigt function

			if (loud)std::cout << "Peak value: " << 2.0 * Vmax_half << "\n"; 

			std::vector<double> dyda(na, 0.0); 
			
			// subtract the peak value since you're looking for point where Voigt = 0.5*Vmax

			Voigt(left, a, &fl, dyda, na); fl -= Vmax_half; // compute the value of the function at the interval endpoints
			Voigt(right, a, &fr, dyda, na); fr -= Vmax_half; // compute the value of the function at the interval endpoints 
			fl = template_funcs::Signum(fl); fr = template_funcs::Signum(fr);

			// test the interval to ensure it contains a root ivttest == -1 => interval has root
			if ( (fl * fr) < 0.0) {

				// the interval contains a root, search can proceed
				if (loud)std::cout << "Initial approximation to the root is " << root << ", fl = "<<fl<<" , fr = "<<fr<<"\n";

				// Compute MAXIT iterations of BisectRoot, stop if the desired tolerance is reached
				count = 0; 
				while (count < MAXIT) {

					// Compute the amount by which the root position must be updated
					dx = 0.5 * (right - left); 

					// Update the position of the root
					root = left + dx; 

					if (loud)std::cout << "Iteration: " << count << ", root = " << root << ", fl = " << fl << " , fr = " << fr << "\n";

					// Test for convergence
					if (fabs(dx) < TOL) {
						if (loud)std::cout << "Bisection has converged to a root within tolerance " << TOL << " after " << count << " iterations\n"; 
						converged = true; 
						break; 
					}
					else {
						// Update the endpoints of the interval containing the root
						Voigt(root, a, &fr, dyda, na); fr -= Vmax_half; fr = template_funcs::Signum(fr); 

						if ( (fl * fr) > 0.0) {
							left = root; 
							fl = fr; 
						}
						else {
							right = root; 
						}
					}

					count++;
				}

				if (converged) {
					if (loud)std::cout << "Root is located at " << root << "\n";
					*HWHM = fabs(root - a[1]);
				}
				else {
					if (loud)std::cout << "Bisection has not converged to a root within tolerance " << TOL << " after " << count << " iterations\n";
					*HWHM = fabs(root - a[1]);
				}
			}
			else {
				// the interval contains no root, search cannot proceed
				std::string reason;
				reason = "Error: testing::Voigt_FWHM()\n";
				reason += "Search interval improperly defined\n";
				throw std::invalid_argument(reason);
			}			
		}
		else {
			std::string reason;
			reason = "Error: testing::Voigt_FWHM()\n";
			if(!c1 || !c2) reason += "Search interval improperly defined\n";
			if (!c3) reason += "Voigt parameters improperly defined\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::Lorentz_Voigt_Fit_Analysis()
{
	// Run the Lorentz and Voigt fitting algorithms for the sample data sets
	// Save the fitting parameters, see if you can fit a model to the Voigt HWHM data
	// R. Sheehan 13 - 12 - 2021

	//std::string data_home = "c:\\users\\robertsheehan\\Research\\Laser_Physics\\Linewidth\\Data\\Sample_Data"; 
	//std::string data_home = "C:\\users\\robertsheehan\\Research\\Laser_Physics\\Linewidth\\Data\\ESA_Spectra_Versus_VOA_Bias\\";
	std::string data_home = "C:\\users\\robertsheehan\\Research\\Laser_Physics\\Linewidth\\Data\\NKT_LCR_DSHI\\";
	
	useful_funcs::set_directory(data_home); 
	useful_funcs::get_directory(); 

	int n_files = 11, npts, n_rows, npars, n_cols, indx_max = 0, ITMAX = 50;
	double spread = 0.15, spctr_max = -500.0, f_max = 0, scale_fac = 1.0e+6, f_start = 60.0, f_end = 100.0, TOL = 1e-3, chisq = 0.0;
	double v_chisq = 0.0, v_nu = 0.0, v_red_chisq = 0.0, v_gof = 0.0; 
	double l_chisq = 0.0, l_nu = 0.0, l_red_chisq = 0.0, l_gof = 0.0; 

	//std::string Vvals[11] = {"000", "100", "200", "300", "325", "350", "360", "365", "370", "375", "400"};
	//double Vvolts[11] = { 0.0, 1.0, 2.0, 3.0, 3.25, 3.5, 3.6, 3.65, 3.7, 3.75, 4.0 }; 

	std::string Vvals[37]; 
	double Vvolts[37]; 

	double dfb = 80, f0 = 80; 
	for (int i = 0; i < 37; i++) {
		Vvolts[i] = f0; Vvals[i] = template_funcs::toString(f0); 
		f0 += dfb; 
	}
	
	std::string file_tmplt; 
	std::string file_out; 
	std::string file_out_gof;
	std::string file_results; 
	
	// Declare the necessary arrays
	std::vector<double> xdata;
	std::vector<double> ydata;
	std::vector<double> sigdata;
	std::vector<std::vector<double>> the_data;

	file_out = "Fitted_Parameter_Values.txt"; 
	//std::ofstream write(file_out, std::ios_base::out , std::ios_base::app);
	std::fstream write(file_out, std::ios_base::out|std::ios_base::trunc);
	if (write.is_open()) {
		write << "Filename , Actual Peak / nW , Voigt h / nW , Voigt f_{0} / MHz , Voigt gamma / MHz , Voigt sigma / MHz , Voigt delta / MHz , Voigt peak / nW , Lorentz h / nW , Lorentz f_{0} / MHz , Lorentz gamma / MHz , Lorentz peak / nW\n"; 
	}

	file_out_gof = "Fitted_Parameter_GOF.txt"; 
	std::fstream write_gof(file_out_gof, std::ios_base::out | std::ios_base::trunc);
	if (write_gof.is_open()) {
		write_gof << "Filename, Voigt chi^{2}, Voigt nu, Voigt chi^{2}/nu, Voigt gof, Lorentz chi^{2}, Lorentz nu, Lorentz chi^{2}/nu, Lorentz gof, chi^{2} ratio\n";
	}

	bool ReScaleFrq = true;
	std::string unitsFrq = ReScaleFrq ? " kHz" : " MHz";
	std::string unitsPow = " pW";

	for (int i = 1; i <=37; i++) {
		//file_tmplt = "Smpl_LLM_" + template_funcs::toString(i) + dottxt; // name of the file with data to be fitted
		//file_tmplt = "LLM_Spctrm_I_" + template_funcs::toString(i) + dottxt; // name of the file with data to be fitted
		//file_tmplt = "JDSU_DFB_T_20_I_50_V_" + Vvals[i-1] + dottxt; 
		file_tmplt = "NKT_I_100_Vb_30_RBW_05_fb_" + Vvals[i-1] + dottxt; 

		vecut::read_into_matrix(file_tmplt, the_data, n_rows, n_cols, true); // read the data from the file into memory

		// format the data from dBm scale to linear scale prior to fitting
		for (int j = 0; j < n_rows; j++) {
			if (ReScaleFrq) {
				xdata.push_back(1000 * (the_data[j][0] - Vvolts[i-1]) ); // rescale frq to zero in units of kHz
			}
			else {
				xdata.push_back(the_data[j][0]); // rescale frq to zero in units of kHz			
			}
			ydata.push_back(scale_fac * convert_dBm_to_mW(the_data[j][1])); // convert the spectral data from dBm to mW scale
		}

		npts = static_cast<int>(xdata.size()); // number of data points in the fit

		// locate the location of the maximum of the data set
		indx_max = 0; spctr_max = -500.0, f_max = 0;
		for (int j = 0; j < npts; j++) {
			if (ydata[j] > spctr_max) {
				spctr_max = ydata[j];
				indx_max = j;
			}

			if (fabs(ydata[j]) > 0.0) {
				sigdata.push_back(spread * ydata[j]);
			}
			else {
				sigdata.push_back(spread); // sigdata cannot have zero values
			}
		}

		f_max = xdata[indx_max]; // find location of peak value

		std::cout << "\nMax value in data set: " << spctr_max << unitsPow << "\n";
		std::cout << "Corresponding Frequency: " << f_max << unitsFrq << "\n\n";

		// Voigt Fit
		// Declare the necessary arrays
		npars = 4; // For Voigt fit 4 parameters are required
		std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
		std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

		// define the initial guesses to the parameters to be determined
		std::vector<double> a_voigt(npars, 0.0);
		std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

		// initial guesses for the parameters	switch to tell code whether or not they're to be fitted
		a_voigt[0] = ReScaleFrq ? spctr_max : 10 * spctr_max;			ia[0] = ReScaleFrq ? 0 : 1;
		a_voigt[1] = f_max;				ia[1] = 0;
		a_voigt[2] = 3.32;				ia[2] = 1;
		a_voigt[3] = 3.32;				ia[3] = 1;

		// run the fitting algorithm
		v_chisq = 0.0; v_nu = 0.0; v_red_chisq = 0.0; v_gof = 0.0; // reset values to 0.0 between calls
		fit::non_lin_fit(xdata, ydata, sigdata, npts, a_voigt, ia, npars, covar, alpha, &v_chisq, Voigt, ITMAX, TOL, true);

		fit::compute_gof(npts, npars, &v_chisq, &v_nu, &v_red_chisq, &v_gof, ia);

		// compute the HWHM using bisection
		double xlow = a_voigt[1], xhigh = xdata[npts - 1], HWHM = 0.0, voigt_peak;

		voigt_peak = a_voigt[0] * exp(template_funcs::DSQR(a_voigt[2] / a_voigt[3])) * probability::erffc(a_voigt[2] / a_voigt[3]); 

		Voigt_HWHM(xlow, xhigh, a_voigt, npars, &HWHM);

		std::cout << "Voigt Fit\nFitted centre freq: " << a_voigt[1] << unitsFrq << "\n";
		std::cout << "Computed peak val: " << voigt_peak << unitsPow << "\n";
		std::cout << "Lorentz HWHM: " << a_voigt[2] << unitsFrq << "\n";
		std::cout << "Gauss HWHM: " << sqrt(2.0 * log(2.0)) * a_voigt[3] << unitsFrq << "\n";
		std::cout << "Voigt HWHM: " << HWHM << unitsFrq << "\n\n";

		// compute the residuals for the fit 
		// Voigt values stored in voigt_data[3], residuals stored in voigt_data[4]
		std::vector<std::vector<double>> voigt_data;
		fit::residuals(xdata, ydata, sigdata, npts, a_voigt, npars, Voigt, voigt_data);

		// clear those vectors between calls to the fit functions
		covar.clear(); alpha.clear(); ia.clear();

		// Lorentz Fit
		// Declare the necessary arrays
		npars = 3; // For Voigt fit 3 parameters are required
		covar = vecut::array_2D(npars, npars);
		alpha = vecut::array_2D(npars, npars);

		// define the initial guesses to the parameters to be determined
		std::vector<double> a_lorentz(npars, 0.0);
		ia.resize(npars, 1); // tell the algorithm that you want to locate all parameters 

		ia[1] = 0;
		a_lorentz[0] = 10 * spctr_max; // must start looking for amplitude parameter 
		a_lorentz[1] = f_max; // this is fixed
		a_lorentz[2] = 1.5; // initial guesses for the parameters

		// run the fitting algorithm
		l_chisq = 0.0; l_nu = 0.0; l_red_chisq = 0.0; l_gof = 0.0; // reset values to 0.0 between calls
		fit::non_lin_fit(xdata, ydata, sigdata, npts, a_lorentz, ia, npars, covar, alpha, &l_chisq, Lorentzian, ITMAX, TOL, true);

		fit::compute_gof(npts, npars, &l_chisq, &l_nu, &l_red_chisq, &l_gof, ia);

		double lorentz_peak = a_lorentz[0] / a_lorentz[2]; 

		std::cout << "Lorentz Fit\nFitted centre freq: " << a_lorentz[1] << unitsFrq << "\n";
		std::cout << "Computed peak val: " << lorentz_peak << unitsPow << "\n";
		std::cout << "Computed HWHM: " << a_lorentz[2] << unitsFrq << "\n";

		// compute the residuals for the fit 
		// Voigt values stored in data[3], residuals stored in data[4]
		std::vector<std::vector<double>> lor_data;
		fit::residuals(xdata, ydata, sigdata, npts, a_lorentz, npars, Lorentzian, lor_data);

		// clear those vectors between calls to the fit functions
		covar.clear(); alpha.clear(); ia.clear();

		// Output fit Parameter Results
		// Output data in the form Filename , Actual Peak / uW , Voigt h / uW , Voigt f_{0} / MHz , Voigt gamma / MHz , Voigt sigma / MHz , Voigt delta / MHz , Voigt peak / uW , Lorentz h / uW , Lorentz f_{0} / MHz , Lorentz gamma / MHz , Lorentz peak / uW
		if (write.is_open()) {
			write << file_tmplt << " , " << std::setprecision(10) << spctr_max << " , "; 
			for (int i = 0; i < 4; i++) {
				write << std::setprecision(10) << a_voigt[i] << " , ";
			}
			write << std::setprecision(10)<< HWHM <<" , " << voigt_peak << " , ";
			for (int i = 0; i < 3; i++) {
				write << std::setprecision(10) << a_lorentz[i] << " , ";
			}
			write << std::setprecision(10) << lorentz_peak << "\n";
		}

		// Output fit gof results
		// Output data in the form Filename, Voigt chi^{2}, Voigt nu, Voigt chi^{2}/nu, Voigt gof, Lorentz chi^{2}, Lorentz nu, Lorentz chi^{2}/nu, Lorentz gof, chi-sq ratio
		if (write_gof.is_open()) {
			write_gof << file_tmplt << " , "; 
			write_gof << std::setprecision(5) << v_chisq << " , " << v_nu << " , " << v_red_chisq << " , " << v_gof << " , " <<l_chisq << " , " << l_nu << " , " << l_red_chisq << " , " << l_gof << " , " << v_chisq / l_chisq << "\n";
		}

		// Output fit function results
		// Data will be output in the form
		// row 1: freq fit data, row 2: LLM spctrl data, row 3: Voigt f vals, row 4: Voigt resids, row 5: Lorentz f vals, row 6: Lorentz f resids
		useful_funcs::remove_substring(file_tmplt, dottxt); 
		file_results = file_tmplt + "_fit_results.txt";
		std::ofstream write_res(file_results, std::ios_base::out, std::ios_base::trunc);
		if (write_res.is_open()) {

			// row 1: freq fit data
			for (int i = 0; i < npts; i++) {
				if (i == npts - 1) {
					write_res << std::setprecision(10) << voigt_data[0][i] << "\n";
				}
				else {
					write_res << std::setprecision(10) << voigt_data[0][i] << " , ";
				}
			}
			
			// row 2: LLM spctrl data
			for (int i = 0; i < npts; i++) {
				if (i == npts - 1) {
					write_res << std::setprecision(10) << voigt_data[1][i] << "\n";
				}
				else {
					write_res << std::setprecision(10) << voigt_data[1][i] << " , ";
				}
			}
			
			// row 3: Voigt f vals
			for (int i = 0; i < npts; i++) {
				if (i == npts - 1) {
					write_res << std::setprecision(10) << voigt_data[3][i] << "\n";
				}
				else {
					write_res << std::setprecision(10) << voigt_data[3][i] << " , ";
				}
			}
			
			// row 4: Voigt resids
			for (int i = 0; i < npts; i++) {
				if (i == npts - 1) {
					write_res << std::setprecision(10) << voigt_data[4][i] << "\n";
				}
				else {
					write_res << std::setprecision(10) << voigt_data[4][i] << " , ";
				}
			}
			
			// row 5: Lorentz f vals
			for (int i = 0; i < npts; i++) {
				if (i == npts - 1) {
					write_res << std::setprecision(10) << lor_data[3][i] << "\n";
				}
				else {
					write_res << std::setprecision(10) << lor_data[3][i] << " , ";
				}
			}
			
			// row 6: Lorentz f resids
			for (int i = 0; i < npts; i++) {
				if (i == npts - 1) {
					write_res << std::setprecision(10) << lor_data[4][i] << "\n";
				}
				else {
					write_res << std::setprecision(10) << lor_data[4][i] << " , ";
				}
			}
			
		}

		write_res.close(); 

		// clear those vectors between file IO
		the_data.clear(); xdata.clear(); ydata.clear(); sigdata.clear(); voigt_data.clear(); lor_data.clear(); 
	}

	write.close(); write_gof.close();
}

void testing::Ring_Down(double x, std::vector<double>& a, double* y, std::vector<double>& dyda, int& na)
{
	// Definition of the model used to approximate a gas ring down process
	// a stores Ring Down parameters a = { A, tau, B}
	// a[0] = A, a[1] = tau, a[2] = B
	// A is an amplitude fitting factor usually expressed in units of mV
	// tau is the Ring Down characterisation time usually expressed in units of us
	// B is an offset from zero parameter usually expressed in units of mV
	// Ring_Down value is given by *y
	// dyda is array that stores value of derivative of Ring_Down function wrt each parameter in a
	// Dimensions of the arrays are a[0..na-1], dyda[0..na-1]
	// na is no. parameters
	// R. Sheehan 16 - 3 - 2021

	try {
		double arg = x / a[1]; // ( t / tau )
		double exp_arg = exp(-1.0 * arg); // exp( -( t / tau ) )
		*y = ( a[0] * exp_arg ) + a[2]; // R = A exp( -(t / tau) ) + B
		dyda[1] = (a[0]/a[1])*arg*exp_arg; // \partial R / \partial tau
		dyda[0] = exp_arg; // \partial R / \partial A
		dyda[2] = 1.0; // \partial T / \partial B
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::Ring_Down_data_fit()
{
	// Apply LM method to measured Ring-Down voltage data
	// R. Sheehan 16 - 3 - 2022

	double A = 3.5; // Amplitude Parameter (mV)
	double tau = 55; // Ring-Down Time Characterisation Parameter ( us )
	double B = 5; // Offset Parameter (mV)
	double t0 = 50; // time after initial decay ( us )
	double y = 0.0; //computed voltage value

	int npars = 3;
	std::vector<double> a(npars, 0.0);
	std::vector<double> dyda(npars, 0.0);

	// Initial guesses for the parameters
	a[0] = A; a[1] = tau; a[2] = B;

	testing::Ring_Down(t0, a, &y, dyda, npars);

	std::cout << "Ring-Down Fit Test Evaluation\n";
	std::cout << "Ring-Down time value: " << t0 << " us\n";
	std::cout << "Ring-Down voltage value: " << y << " V\n";
	std::cout << "Der wrt A: " << dyda[0] << "\n";
	std::cout << "Der wrt tau: " << dyda[1] << "\n";
	std::cout << "Der wrt B: " << dyda[2] << "\n\n";
	std::cout << "Rign-Down Equation Fit\n";

	// Generate data to use in the fit process
	int npts;
	long idum = (-1011);
	double spread, xlow, xhigh, deltax, xpos, yval;

	// Generate initial data set
	xlow = 0; xhigh = 300.0; deltax = 1.0;
	npts = (int)(1.0 + ((xhigh - xlow) / deltax));

	std::vector<double> xdata(npts, 0.0);
	std::vector<double> ydata(npts, 0.0);
	std::vector<double> sigdata(npts, 0.0);

	spread = 0.1; // variance of the noise being added to the signal
	xpos = xlow;
	for (int i = 0; i < npts; i++) {

		Ring_Down(xpos, a, &yval, dyda, npars); // evaluate the diode-voltage function

		xdata[i] = xpos;

		yval *= rng::gasdev1(&idum, 1.0, template_funcs::DSQR(spread)); // add noise to the signal value

		ydata[i] = yval;

		sigdata[i] = fabs(yval) > 0.0 ? spread * yval : spread; // sigdata cannot have zero values

		xpos += deltax;
	}

	// Perform the best it search for the data set
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	//ia[1] = 0; // search for params 0 and 2, fix param 1 value
	a_guess[0] = 7.0; a_guess[1] = 50.0; a_guess[2] = 0.0; // initial guesses for the parameters

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, Ring_Down, ITMAX, TOL, true);

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, Ring_Down, data);

	// output the residuals 
	std::string thefile = "Ring_Down_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	a.clear();
}

void testing::Ring_Down_data_fit_test()
{
	// Apply LM method to measured ring down data
	// R. Sheehan 16 - 3 - 2022

	// Read in the measured spectral data
	//std::string filename = "LC_LT372_PROC_22_03_2022_12_02_M0.txt";
	std::string filename = "LC_LT372_PROC_22_03_2022_13_02_M0.txt";
	//std::string filename = "LC_LT372_PROC_22_03_2022_13_52_M0.txt";

	int npts, n_rows, npars = 3, n_cols;
	double spread = 0.15;

	time_t rawtime;
	time(&rawtime);
	long idum = (-static_cast<long>(rawtime)); // randomised seed for rng

	std::vector<std::vector<double>> the_data;

	vecut::read_into_matrix(filename, the_data, n_rows, n_cols, '#', true);

	// Change the units of the imported data time / us, voltage / mV
	for (int i = 0; i < n_rows; i++) {
		the_data[i][0] *= 1.0e+6; 
		the_data[i][1] *= 1.0e+3; 
	}

	std::cout << "The Data\n"; 
	vecut::print_mat(the_data);

	// store data locally
	std::vector<double> xdata;
	std::vector<double> ydata;
	std::vector<double> sigdata;

	npts = 1024;

	for (int i = 0; i < npts; i++) {
		xdata.push_back(the_data[i][0]); 
		ydata.push_back(the_data[i][1]); 
		if (fabs(the_data[i][1]) > 0.0) {
			sigdata.push_back(spread * the_data[i][1]);
		}
		else {
			sigdata.push_back(spread); // sigdata cannot have zero values
		}
	}

	//npts = static_cast<int>(xdata.size()); 

	// Perform the best it search for the data set
	int ITMAX = 10;

	double TOL = 0.001;
	double chisq = 0.0;

	// Declare the necessary arrays
	std::vector<std::vector<double>> covar = vecut::array_2D(npars, npars);
	std::vector<std::vector<double>> alpha = vecut::array_2D(npars, npars);

	// define the initial guesses to the parameters to be determined
	std::vector<double> a_guess(npars, 0.0);
	std::vector<int> ia(npars, 1); // tell the algorithm that you want to locate all parameters 

	a_guess[0] = 7.0; a_guess[1] = 50.0; a_guess[2] = 0.2; // initial guesses for the parameters

	// run the fitting algorithm
	fit::non_lin_fit(xdata, ydata, sigdata, npts, a_guess, ia, npars, covar, alpha, &chisq, Ring_Down, ITMAX, TOL, true);

	// compute the residuals for the fit
	std::vector<std::vector<double>> data;
	fit::residuals(xdata, ydata, sigdata, npts, a_guess, npars, Ring_Down, data);

	// output the residuals 
	std::string thefile = "Ring_Down_non_lin_fit.txt";

	int nrows = 5;
	vecut::write_into_file(thefile, data, nrows, npts);

	// take a look at the goodness of fit statistics
	double chisqr = 0.0, rsqr = 0.0, dof = static_cast<int>(npts - npars), gof = 0.0;
	fit::goodness_of_fit(xdata, ydata, data[4], npts, a_guess, npars, Ring_Down, &chisqr, &dof, &rsqr, &gof);

	xdata.clear(); ydata.clear(); sigdata.clear();
	a_guess.clear(); ia.clear(); covar.clear(); alpha.clear();
	data.clear(); the_data.clear();
}