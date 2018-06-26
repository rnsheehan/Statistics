#ifndef ATTACH_H
#include "Attach.h"
#endif

// The aim of this project is to implement the Statistical routines described in NRinC, Ch14

// The aim is to create a namespace holding all the statistical tests that are described. 
// Methods described to date include
// moment, avevar, ttest, tutest, tptest, ftest

// Also required is the implementation of a namespace that computes values of the functions
// needed to compute probabilities, some of these are already in existence, some are not. 
// Create a namespace to hold the required probability functions

// Use std::vector and exception handling where appropriate

// R. Sheehan 5 - 9 - 2017

void Gamm_Test(); 

void sort_test();
void sort2_test(); 
//void sort2(std::vector<double> &arr1, std::vector<double> &arr2); 
//void test_func(std::vector<double> &name); 

void sand_box(); 

int main()
{
	//testing::compute_moments(); 

	//testing::t_test_diff_mean_same_var(); 

	//testing::t_test_diff_mean_diff_var(); 

	//testing::f_test_diff_mean_diff_var(); 

	//testing::chsone_test(); 

	//testing::ksone_test(); 

	//testing::kstwo_test();

	//sort2_test(); 

	/*std::string Head; 

	std::vector<std::string> Mths; 

	std::vector<std::string> Acc; 

	std::vector<std::string> Cities;

	std::vector<std::string> States;

	std::vector<int> State_Averages; 

	std::vector<double> State_Latitudes; 

	std::vector< std::vector<int> > Num;
	
	testing::read_tab1_data(Head, Mths, Acc, Num); 

	testing::read_tab2_data(Head, Mths, Cities, States, State_Averages, State_Latitudes, Num); */

	//std::string test = "Testing the string\n"; 
	//if(!(test.find("xyz") != std::string::npos)) std::cout<<test; 

	//testing::con_tab_test(); 

	//testing::pearson_test(); 

	//testing::spearman_test(); 

	//testing::rank_order_test(); 

	//testing::kendall_test_1(); 

	//testing::kendall_test_2();

	//testing::lin_fit_test(); 

	//testing::gaussj_test(); 

	testing::lst_sqr_test(); 
	
	std::cout<<"Press enter to close console\n";
	std::cin.get(); 

	return 0; 
}

void Gamm_Test()
{
	// output some values of the gamma related functions

	double x = 2.0; 
	double a = 0.25; 

	std::cout<<probability::gammln(x)<<"\n";
	std::cout<<probability::gammp(a,x)<<"\n";
	std::cout<<probability::gammq(a,x)<<"\n";

	double galn, gaser, gacf;

	probability::gser(&gaser, a, x, &galn); 

	std::cout<<gaser<<"\n"; 
	std::cout<<galn<<"\n"; 

	probability::gcf(&gacf, a, x, &galn); 

	std::cout<<gacf<<"\n"; 
	std::cout<<galn<<"\n"; 
}

void sort_test()
{
	// Check the operation of the std::sort algorithm
	// http://www.cplusplus.com/reference/algorithm/sort/
	// R. Sheehan 8 - 9 - 2017

	std::vector<double> data; 

	data.push_back(-9); data.push_back(7); data.push_back(22); data.push_back(-3);
	data.push_back(12); data.push_back(77); data.push_back(-102); data.push_back(23); 
	data.push_back(6); data.push_back(1); data.push_back(11); data.push_back(-5); 
	data.push_back(-17); data.push_back(3); data.push_back(17); data.push_back(-5); 

	std::cout<<"The data is\n";
	for(int i=0; i<(int)(data.size()); i++)std::cout<<data[i]<<"\n"; 

	std::sort(data.begin(), data.end()); // sorts in ascending order by default

	std::cout<<"\nThe sorted data is\n";
	//for(int i=0; i<(int)(data.size()); i++)std::cout<<data[i]<<"\n";
	
	std::vector<double>::iterator it = data.begin(); 
	for(it; it != data.end(); it++)std::cout<<*it<<"\n"; 

	data.clear(); 
}

void sort2_test()
{
	// In NRinC a function sort2 is implemented which sorts a column of data and performs the same re-arrangements on a second column
	// The algorithm is used in the calculation of the Spearman correlation coefficient
	// I'd like to use std::sort to do the sorting if possible. 
	// R. Sheehan 9 - 1 - 2018

	// Make two columns of data
	std::vector<double> data1, data2; 

	data1.push_back(-9); data1.push_back(7); data1.push_back(22); data1.push_back(-3);
	data1.push_back(12); data1.push_back(77); data1.push_back(-102); data1.push_back(23); 
	data1.push_back(6); data1.push_back(1); data1.push_back(11); data1.push_back(-5); 
	data1.push_back(-17); data1.push_back(3); data1.push_back(17); data1.push_back(-5); 

	data2.push_back(-6); data2.push_back(-1); data2.push_back(-11); data2.push_back(-5); 
	data2.push_back(9); data2.push_back(-7); data2.push_back(-22); data2.push_back(3);
	data2.push_back(17); data2.push_back(-3); data2.push_back(17); data2.push_back(5);
	data2.push_back(-12); data2.push_back(-77); data2.push_back(102); data2.push_back(-23); 

	// Join them together to form a single array
	// use the pair data type rather than std::vector<std::vector<double>> Data; 
	std::vector< std::pair< double, double > > Data; 

	std::vector<double>::iterator it1 = data1.begin(); 
	std::vector<double>::iterator it2 = data2.begin(); 
	for(it1; it1 != data1.end(); it1++){
		Data.push_back( std::make_pair(*it1, *it2) ); 
		it2++; 
	}

	std::cout<<"The data is\n";
	for(int i=0; i<(int)(Data.size()); i++){
		std::cout<<Data[i].first<<" , "<<Data[i].second<<"\n"; 
	}

	std::sort(Data.begin(), Data.end()); // looks like it will sort according to the first element of the pair

	std::cout<<"\nThe sorted data is\n";
	for(int i=0; i<(int)(Data.size()); i++){
		std::cout<<Data[i].first<<" , "<<Data[i].second<<"\n"; 
	}

	// Unpack the data into the original vectors
	std::cout<<"\nThe sorted data is\n";
	for(int i=0; i<(int)(Data.size()); i++){
		data1[i] = Data[i].first; 
		data2[i] = Data[i].second; 
		std::cout<<data1[i]<<" , "<<data2[i]<<"\n"; 
	}

	// try the sort2 function
	//sort2(data2, data1); 
	template_funcs::sort2(data2, data1); 

	std::cout<<"\nThe sorted data is\n";
	for(int i=0; i<(int)(data2.size()); i++){
		std::cout<<data1[i]<<" , "<<data2[i]<<"\n"; 
	}


	data1.clear(); data2.clear(); Data.clear(); 
}

void sand_box()
{
	/*bool c1 = false; 
	bool c2 = true; 

	if (c1) std::cout << "c1 is true\n"; 
	if (!c1) std::cout << "c1 is false\n"; 
	if (c2) std::cout << "c2 is true\n";
	if (!c2) std::cout << "c2 is false\n";*/

	int rows = 9; 
	int cols = 7; 
	std::vector<std::vector<double>> X; 

	X = lin_alg::array_2D(rows, cols); 

	std::pair<int, int> test = lin_alg::array_2D_size(X); 

	X.clear(); 
}

//void sort2(std::vector<double> &arr1, std::vector<double> &arr2)
//{
//	// Sort arr1 into ascending order using while making the corresponding rearrangement of the arr2
//	// This is an implementation of the NRinC::sort2 using std::sort
//	// arr1 and arr2 must have the same number of elements
//	// R. Sheehan 10 - 1 - 2018
//	
//	try{
//		bool c1 = arr1.size() == arr2.size() ? true : false; // arr1 and arr2 must have the same dimensions!
//		bool c2 = arr1.size() > 0 ? true : false; 
//		bool c3 = arr2.size() > 0 ? true : false; 
//		bool c4 = (c1 && c2 && c3) ? true : false; 
//		if(c4){
//			// Join the data sets together to form a single container
//			// use the pair data type rather than std::vector<std::vector<double>> Data; 
//			std::vector< std::pair< double, double > > Data; 
//
//			// Store the data from arr1 and arr2 in Data
//			// arr1 and arr2 must have the same dimensions!
//			// arr1 corresponds to the first column, arr2 to the second
//			std::vector<double>::iterator it1 = arr1.begin(); 
//			std::vector<double>::iterator it2 = arr2.begin(); 
//			for(it1; it1 != arr1.end(); it1++){
//				Data.push_back( std::make_pair(*it1, *it2) ); 
//				it2++; 
//			}
//
//			// sort the data using std::sort
//			// this will sort according to the data in the first element of the pair
//			// the same rearrangement will be applied to the data in the second element
//			std::sort(Data.begin(), Data.end()); 
//
//			// place the sorted data back into the original containers		
//			/*for(int i=0; i<(int)(Data.size()); i++){
//				arr1[i] = Data[i].first; arr2[i] = Data[i].second; 
//			}*/		
//			std::vector<std::pair<double,double>>::iterator pit = Data.begin(); 
//			int i=0; 
//			for(pit; pit != Data.end(); pit++){
//				arr1[i] = pit->first; arr2[i] = pit->second; 
//				i++; 
//			}
//
//			Data.clear(); 
//		}
//		else{
//			std::string reason; 
//			reason = "Error::sort2\n"; 
//			if(!c1) reason+="arr1 and arr2 have different sizes\n"; 
//			if(!c2) reason+="arr1 has no elements\n";
//			if(!c3) reason+="arr1 has no elements\n";
//			throw std::invalid_argument(reason); 
//		}
//	}
//	catch(std::invalid_argument &e){
//		std::cerr<<e.what(); 
//	}
//}

