# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 14:15:40 2018

@author: am17010
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import PowerNorm
import scipy
from scipy import integrate

k = 6.3e6                                                                       # declaring the parameters that will remain constant
x1 = -0.005
x2 = 0.005

def function(xPrime, k, z, x):                                                  # defining the function that contains the Fresnel integral
                                                                                # k is the wavenumber, z is the aperture-screen distance, x and y values are screen coordinates, and xPrime and yPrime values are aperture coordinates
    return math.e ** (1j * k / (2 * z) * (x - xPrime) ** 2)



def Simpson(x1Prime, x2Prime, k, z, x, N):                                      # defining the function that uses Simpson's rule to approximate the Fresnel integral
    
    deltax = (x2Prime - x1Prime) / N
    sumOdd = 0                                                                  # declaring the odd and even sums 
    sumEven = 0
    
    x1DoublePrime = x1Prime + deltax                                            # the double primes act as dummy variables so the originals, x1Prime and x2Prime are not overwritten
    x2DoublePrime = x1Prime + 2 * deltax
    
    for i in range(1, N // 2 + 1):
        
        sumOdd += 4 * function(x1DoublePrime, k, z, x)
        x1DoublePrime += 2 * deltax
    
    for i in range(1, N // 2):
        
        sumEven += 2 * function(x2DoublePrime, k, z, x)
        x2DoublePrime += 2 * deltax
    
    total = (k / (2 * math.pi * z)) * (sumOdd + sumEven +                       # summing the separate parts of the equation to get the full approximation
            function(x1Prime, k, z, x) + function(x2Prime, k, z, x)) * (deltax 
                    / 3)
    total = abs(total) ** 2
    
    return total



def valueCheckPosFloat(value):                                                  # this function checks that the user has input a positive float when they're meant to
    
    while value.isalpha() or float(value) <= 0:
        
        print("You require a positive float!")
        value = input("Please enter a positive float: ")
    
    value = float(value)
    
    return value
    


def valueCheckN(value):                                                         # this function checks that the user has input a positive, even integer when they're meant to
    
    while value.isalpha() or float(value) <= 0 or float(value) % 2 != 0:
        
        print("N must be a positive, even integer!")
        value = input("Enter the number of intervals (even integer): ")
   
    value = int(value)
    
    return value



def valueCheckGraphs(value):                                                    # this function checks that the user has input a positive integer greater than 1 when they're meant to
    
    while value.isalpha() or float(value) <= 1 or float(value) % 1 != 0:
        
        print("The number of graphs must be a positive integer greater than 1!")
        value = input("Enter the number of graphs you want: ")
    
    value = int(value)
    
    return value



def integral(xPrimemax, yPrimemax, x, y, z, k):                                 # this function uses scipy to calculate the integral at the x and y values on the screen
                                                                                
    xvals = []
    yvals = []
        
    xPrimevals = np.linspace(x1Prime, x2Prime, N)
    yPrimevals = np.linspace(y1Prime, y2Prime, N)
            
    for i in xPrimevals:                                                        # calculating the function inside the integral at each value on the slit
        
        xvals.append(function(i, k, z, x))
        
    for i in yPrimevals:
        
        yvals.append(function(i, k, z, y))
                
    xScreenval = scipy.integrate.simps(xvals, xPrimevals)                       # using scipi to find the integral
    yScreenval = scipy.integrate.simps(yvals, yPrimevals)        
    intensity = (abs(xScreenval * yScreenval)) ** 2     
    
    return intensity
       


def image(x2Prime, y2Prime, screenWidth, z, k, N):                              # this function produces the image by creating an N x N array and applies the integral to each coordinate
    
    image = np.zeros((N, N))
    xycoord = np.linspace(-(screenWidth / 2), (screenWidth / 2), N)             # defines where to start and stop the image (+/- half the screen width from the centre
            
    for i in range(N):
        for j in range(N):
                    
            image[i, j] = integral(xPrimemax, yPrimemax, xycoord[i], 
                 xycoord[j], z, k)
                    
    image = (k / (2 * math.pi * z)) * image
    
    return image
    


MyInput = "0"

while MyInput != "q":                                                           # menu system for ease of use
    
    print("\n")
    print("Option a creates a 1D intensity graph with variables input by the user.")
    print("\n")
    print("Option b is used to compare 1D intensity graphs with varying distances.")
    print("\n")
    print("Option c creates a 2D intensity heatmap with the ability to vary")
    print("the slit width in the y direction as well as x.")
    print("\n")
    
    MyInput = input("Enter a choice, 'a', 'b', 'c', or 'q' to quit: ")
    print("You entered the choice: " + MyInput)
    
    
    
    if MyInput == "a":                                                          # part a explained in the code
        
        print("\n")                                                            
        print("This part asks the user (you) for the slit-screen distance,") 
        print("the slit width, and the number of intervals. It then produces") 
        print("a graph with these variables.")
        
        N = input("Enter the number of intervals (even integer): ")
        N = valueCheckN(N)

        z = input("Enter the distance to the screen (~0.02): ")
        z = valueCheckPosFloat(z)
   
        xPrimemax = input("Enter the slit width (~2e-5): ")
        xPrimemax = valueCheckPosFloat(xPrimemax)
        x2Prime = xPrimemax / 2
        x1Prime = - (xPrimemax) / 2
    
        i = x1
        yvals = []                                                              # creating empty lists to store the x and y values that will be used in the graph
        xvals = []
        
        while i <= x2:
            
            yvals.append(Simpson(x1Prime, x2Prime, k, z, i, N))
            xvals.append(i)
            i += (x2 - x1) / N

        plt.plot(xvals, yvals)
        plt.xlabel("Screen coordinate (m)")
        plt.ylabel("Intensity")
        plt.show()
        
    
    
    elif MyInput == "b":                                                        # part b also explained in the code
        
        print("\n")                                                            
        print("This part asks you for the number of graphs you want plotted")
        print("and maximum and minimum values for the slit-screen distance")
        print("and the slit width. These graphs are then plotted going from")
        print("the minimum to maximum of each value.")
        print("(If you want to keep the distance or width constant, simply")
        print("use the same value for the maximum and minimum).")
        
        noGraphs = input("Enter the number of graphs you want: ")
        noGraphs = valueCheckGraphs(noGraphs)
        
        N = input("Enter the number of intervals (even integer): ")
        N = valueCheckN(N)
        
        zRange1 = input("Enter the minimum value for the screen distance (~0.02): ")
        zRange1 = valueCheckPosFloat(zRange1)
        
        zRange2 = input("Enter the maximum value for the screen distance (~0.02): ")
        zRange2 = valueCheckPosFloat(zRange2)
        
        slitWidth1max = input("Enter the minimum value of the slit width (~2e-5): ")
        slitWidth1max = valueCheckPosFloat(slitWidth1max)
        slitWidth1 = - (slitWidth1max / 2)
        slitWidth2 = slitWidth1max / 2
        
        slitWidth3max = input("Enter the maximum value of the slit width (~2e-5): ")
        slitWidth3max = valueCheckPosFloat(slitWidth3max)
        slitWidth3 = - (slitWidth3max / 2)
        slitWidth4 = slitWidth3max / 2
        
        slitWidthvals2 = []                                                     # declaring the empty lists that will contain the values of z and the slit width
        slitWidthvals1 = []  
        zvals = []
        
        
        zRange = (zRange2 - zRange1) / (noGraphs - 1)
        
        for i in range(noGraphs + 1):
            
            zvals.append(zRange1 + i * zRange)
            
        slitWidthRange = ((slitWidth4 - slitWidth3) - (slitWidth2 - 
                          slitWidth1)) / (noGraphs - 1)
        
        
        for i in range(noGraphs + 1):
            
            slitWidthvals2.append(slitWidth2 + (i / 2) * slitWidthRange)
            slitWidthvals1.append(slitWidth1 - (i / 2) * slitWidthRange)
            
        for i in range(noGraphs):                                               # for loop so every graph is printed on the same axes
            
            p = x1
            yvals = []
            xvals = []
        
            while p <= x2:
            
                yvals.append(Simpson(slitWidthvals1[i], slitWidthvals2[i], k, 
                                     zvals[i], p, N))
                xvals.append(p)
                p += (x2 - x1) / N

            plt.plot(xvals, yvals, label = "z = " + str(round(zvals[i], 4)) + 
                     ", slit width = " + str(round((slitWidthvals2[i] - 
                                                    slitWidthvals1[i]), 7)))
            plt.legend(loc = 1 , prop = {'size': 8})
            plt.xlabel("Screen coordinate (m)")
            plt.ylabel("Intensity")
        
        plt.show()
        
    
    
    elif MyInput == "c":                                                        # part c explained in the code and functions defined at the beginning
        
        print("\n")
        print("This part uses the inbuilt scipy Simpson function to display")
        print("a heatmap of the intensity through the use of a 2D array, and")
        print("allows the user to vary the aspect ratio of the slit width.")
        
        N = input("Enter the number of intervals (even integer): ")
        N = valueCheckN(N)
        
        z = input("Enter the distance to the screen (~0.02 for far-field, ~5e-3 for near-field): ")
        z = valueCheckPosFloat(z)
   
        xPrimemax = input("Enter the slit width in the x direction (~2e-5 for far-field, ~4e-4 for near-field): ")
        xPrimemax = valueCheckPosFloat(xPrimemax)
        x2Prime = xPrimemax / 2
        x1Prime = - (xPrimemax) / 2
        
        yPrimemax = input("Enter the slit width in the y direction (~2e-5 for far-field, ~4e-4 for near-field): ")                                                                               
        yPrimemax = valueCheckPosFloat(yPrimemax)                               # creating a yPrime variable to vary the slit width in the y direction
        y2Prime = yPrimemax / 2
        y1Prime = - (yPrimemax) / 2
               
        screenWidth = input("Enter the screen width (~9e-3 for far-field, ~8e-4 for near-field): ")
        screenWidth = valueCheckPosFloat(screenWidth)
        
        graph = image(x2Prime, y2Prime, screenWidth, z, k, N)
           
        plt.imshow(graph, cmap = cm.hot, norm = PowerNorm(0.6),                 # changing the gamma to make the fringes more visible
                   origin = "lower", extent = [-graph.shape[1] *                # changing the coordinates on the graph so they match those expected
                                               (screenWidth / (4 * N)), 
                                               graph.shape[1] * 
                                               (screenWidth / (4 * N)), 
                                               -graph.shape[0] * 
                                               (screenWidth / (4 * N)), 
                                               graph.shape[0] * 
                                               (screenWidth / (4 * N))])             
        plt.xlabel("Screen coordinate (m)")
        plt.ylabel("Screen coordinate (m)")
        plt.show()
        
    elif MyInput == "q":
        
        print("You have chosen to quit.")
    
    else:
        
        print("Please enter a valid input!")
    



    









