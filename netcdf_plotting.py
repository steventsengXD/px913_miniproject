#################    PX913 MINI PROJECT      #################

#This part was written by Charlotte Rogerson
#For the PX913 Mini Project

#This code imports the NetCDF file and plots
#the Ex variable and plots the particle position
#in both x and y.

#Uploading the NetCDF file
import matplotlib.pyplot as plt
import netCDF4 as nc

#Accessing the dataset
file = nc.Dataset('electron.nc4', mode='r', format="NETCDF4") 
#Sanity checks - uncomment if not needed
print(file)

#Defining the variables in the file for the pseudocolour plot
#Only need Ex, Ey - Electric fields
Ex = file.variables['Ex_field']
Ey = file.variables['Ey_field']

#Variables for the scatter plot
#Position
pos = file.variables['Position']

#Sorting out the axis limits
ymin = min(pos[1,:])
ymax = max(pos[1,:])

xmin = min(pos[0,:])
xmax = max(pos[0,:])

#Pseudocolour plot of the variable Ex
#If it is being too slow use instead pcolormesh
#Need a 2d grid
plt.figure(1)

plt.pcolor(Ex, cmap = 'Blues' )
plt.title('Electric Field (Ex)')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

#Scatter graph of particle posion in x vs y
plt.figure(2)

plt.scatter(pos[0,:], pos[1,:])
plt.title('Electron Position')
plt.xlabel('Electron position in X')
plt.ylabel('Electron position in Y')
plt.ylim(ymin,ymax)
plt.xlim(xmin, xmax)
plt.annotate("Starting point", (pos[0,0], pos[1,0]))
plt.show()