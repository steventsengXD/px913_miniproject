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
file = nc.Dataset('electron.nc4', mode='r', format="NETCDF4")
#Sanity checks - uncomment if not needed
# print(file)

#Defining the variables in the file for the pseudocolour plot
#Only need Ex, Ey - Electric fields
Ex = np.array(file.variables['Ex_field'])
Ey = np.array(file.variables['Ey_field'])

#Sanity checks - uncomment if not needed
# print(Ey)
# print(Ex)

#Variables for the scatter plot
#Position
pos = np.array(file.variables['Position'])

#Sanity checks - uncomment if not needed
# print(pos)
#pos[:,0] = x_data, pos[:,1] = y_data,

#This sets the grid size
Nx=Ex.shape[0]
Ny=Ey.shape[0]
xdata = np.linspace(-1, 1, Nx)
ydata = np.linspace(-1, 1, Ny)

#This finds the proper bounds for the position plot
if abs(np.max(pos[1,:])-np.min(pos[1,:])) < 0.01:
    ylb=np.min(pos[1,:])*1.25
    yub=np.max(pos[1,:])*1.25
    xlb=np.min(pos[0,:])*1.15
    xub=np.max(pos[0,:])*1.15
else:
    ylb=np.min(pos[1,:])-0.05
    yub=np.max(pos[1,:])+0.05
    xlb=np.min(pos[0,:])-0.05
    xub=np.max(pos[0,:])+0.05



#Pseudocolour plot of the variable Ex
#If it is being too slow use instead pcolormesh
#Need a 2d grid
fig,ax=plt.subplots(nrows=1,ncols=2,figsize=(13,5))

plt.subplot(1,2,1)
plt.pcolor(xdata, ydata, Ex, cmap = 'Blues') #, shading='auto')
plt.title('Electric Field (Ex)')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.colorbar()

#Scatter graph of particle positioon in x vs y
plt.subplot(1,2,2)
plt.scatter(pos[0,1:], pos[1,1:])
plt.scatter(pos[0,0], pos[1,1], label='Initial Position')
plt.title('Electron Position')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.ylim(ylb,yub)
plt.xlim(xlb,xub)
plt.legend()

plt.tight_layout()
plt.show()




# plt.figure(1)
#
# plt.pcolor(xdata, ydata, Ex, cmap = 'Blues' )
#
# plt.title('Electric Field (Ex)')
# plt.xlabel('X')
# plt.ylabel('Y')
#
# plt.show()
#
#
# #Scatter graph of particle positioon in x vs y
# plt.figure(2)
#
# plt.scatter(pos[:,0], pos[:,1])
# plt.title('Electron Position')
# plt.xlabel('Electron position in X')
# plt.ylabel('Electron position in Y')
#
# plt.show()
