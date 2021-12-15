# Import libraries
# You should try an import the bare minimum of modules
import sys # access system routines
import os
import glob
import re

import math
import scipy
import numpy as np
import matplotlib.pyplot as plt

# add path to our file
sys.path.append('c:/Users/robertsheehan/Programming/Python/Common/')
sys.path.append('c:/Users/robertsheehan/Programming/Python/Plotting/')

import Common
import Plotting

MOD_NAME_STR = "Plots" # use this in exception handling messages
       
def non_lin_fit_plot(filename):
    # make a plot of the data from the non-linear fit tests
    # R. Sheehan 21 - 10 - 2021

    FUNC_NAME = ".non_lin_fit_plot()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
            
        if glob.glob(filename):
            # import the dataset
            hv_data = []; labels = []; marks = []; 
            hv_data_1 = []; labels_1 = []; marks_1 = []; 
            data = np.loadtxt(filename, delimiter = ',')
            hv_data.append([data[0], data[1]]); marks.append(Plotting.labs_pts[0]); labels.append('data'); 
            hv_data.append([data[0], data[3]]); marks.append(Plotting.labs_lins[1]); labels.append('fit'); 
            hv_data_1.append([data[0], data[2]]); marks_1.append(Plotting.labs_pts[1]); labels_1.append('sigdata'); 
            hv_data_1.append([data[0], data[4]]); marks_1.append(Plotting.labs_pts[2]); labels_1.append('residuals'); 
            
            # plot the original data with the fitted function                
            args = Plotting.plot_arg_multiple()
            
            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'X'
            args.y_label = 'Y'
            args.fig_name = filename.replace('.txt','')
            args.plt_title = filename.replace('.txt','')
            
            Plotting.plot_multiple_curves(hv_data, args)
            
            # plot the original data with the fitted function                
            args = Plotting.plot_arg_multiple()
            
            args.loud = True
            args.crv_lab_list = labels_1
            args.mrk_list = marks_1
            args.x_label = 'X'
            args.y_label = 'Y'
            args.fig_name = filename.replace('.txt','') + '_Resid'
            args.plt_title = filename.replace('.txt','') + '_Resid'
            
            #Plotting.plot_multiple_curves(hv_data_1, args)
            
        else:
            ERR_STATEMENT = ERR_STATEMENT + "\nFile: " + filename + " not found"
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def main():
    pass

if __name__ == '__main__':
    main()

    pwd = os.getcwd() # get current working directory
    
    print(pwd)
    
    filename = "Diode_non_lin_fit.txt"
    non_lin_fit_plot(filename)
    
    #filename = "Gauss_non_lin_fit.txt"
    #filename = 'Diode_non_lin_fit.txt'
    
    #filename = 'Lorentzian_non_lin_fit.txt'    
    #non_lin_fit_plot(filename)
    
    #filename = 'Gaussian_non_lin_fit.txt'
    #non_lin_fit_plot(filename)

    #filename = 'Voigt_non_lin_fit.txt'
    #non_lin_fit_plot(filename)
    
