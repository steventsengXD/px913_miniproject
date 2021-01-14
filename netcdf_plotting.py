#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 18:36:11 2021

@author: charlotte
"""

#################    PX913 MINI PROJECT      #################

#This part was written by Charlotte Rogerson
#For the PX913 Mini Project

#This code imports the NetCDF file and plots
#the Ex variable and plots the particle position
#in both x and y.


#Uploading the NetCDF file
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

#Accessing the dataset
file = nc.Dataset('filename.nc4', mode='r', format="NETCDF4") 
#reads in the file

#Defining the variables in the file for the pseudocolour plot
#Only need Ex, Ey - Electric fields
Ex = file.variables['Ex_field']
Ey = file.variables['Ey_field']

#Sanity checks - uncomment if not needed
print(Ey)
print(Ex)

#Variables for the scatter plot
#Position
pos = file.variables['Position']

#Sanity checks - uncomment if not needed
print(pos)

#This sets the grid size - possibly change later to be user defined
xdata = np.linspace(-1, 1, 100)
ydata = np.linspace(-1, 1, 100)

#Sanity checks - uncomment if not needed
print(file)


#Pseudocolour plot of the variable Ex
#If it is being too slow use instead pcolormesh
#Need a 2d grid
plt.figure(1)
plt.title('Electric Field (Ex)')
plt.xlabel('X')
plt.ylabel('Y')

plt.show()


#Scatter graph of particle posion in x vs y
plt.figure(2)
plt.scatter(x,y)

plt.xlabel('Particle position in X')
plt.ylabel('Particle position in Y')

plt.show()