# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 11:01:49 2017

@author: Robert Sheehan

Use this module to generate test random number data based on specific distributions
The numpy module is better developed than anything I'm going to write so it has access to more distributions
"""

import os
import numpy
import matplotlib.pyplot as plt

def write_to_file(file_name, numbers):
    # write a generated set of random numbers to a file
    try:
        if numbers is not None:
            thefile = open(file_name,"w")
            
            if thefile.closed:
                raise Exception
            else:
                for x in numbers:
                    thefile.write("%(num)0.9f\n"%{"num":x})
                    
                thefile.close()
        else:
            raise Exception
    except Exception:
        print("\nError: write_to_file()\n")
    
def plot_histogram(data, the_title):
    # make a basic histogram plot of some dat
    
    try:
        if data is not None:
            n_bins = max( len(data)/100, 10)
            print(n_bins)
            plt.hist(data, bins = n_bins, normed = True)
            plt.title(the_title)
            plt.xlabel("Value")
            plt.ylabel("Frequency")
            plt.show()
            plt.clf()
            plt.cla()
            plt.close()
    except Exception:
        print("\nError: plot_histogram()\n")

def main():
    pass

if __name__ == '__main__':
    main()

    print(os.getcwd())
    
    file1 = "Data_1.txt"
    mean = 0.0
    var = 1.0
    size = 1000
    data1 = numpy.random.normal(mean, var, size)
    
    plot_histogram(data1, file1)
    
    write_to_file(file1, data1)
    
    file1 = "Data_2.txt"
    mean = 1.0
    var = 0.2
    size = 1000
    data1 = numpy.random.normal(mean, var, size)
    
    plot_histogram(data1, file1)
    
    write_to_file(file1, data1)
    
    file1 = "Data_3.txt"
    mean = 0.5
    var = 0.5
    size = 1000
    data1 = numpy.random.normal(mean, var, size)
    
    plot_histogram(data1, file1)
    
    write_to_file(file1, data1)
    
    file1 = "Data_3.txt"
    mean = 1.5
    var = 0.5
    size = 1000
    data1 = numpy.random.normal(mean, var, size)
    
    plot_histogram(data1, file1)
    
    write_to_file(file1, data1)
    
    file1 = "Data_4.txt"
    mean = 0.0
    var = 0.5
    size = 1000
    data1 = numpy.random.normal(mean, var, size)
    
    plot_histogram(data1, file1)
    
    write_to_file(file1, data1)
    
    file1 = "Data_5.txt"
    mean = 0.01
    var = 0.5
    size = 1000
    data1 = numpy.random.normal(mean, var, size)
    
    plot_histogram(data1, file1)
    
    write_to_file(file1, data1)
    
    file1 = "Data_6.txt"
    mean = 1.0
    var = 1.0
    size = 1000
    data1 = numpy.random.normal(mean, var, size)
    
    plot_histogram(data1, file1)
    
    write_to_file(file1, data1)