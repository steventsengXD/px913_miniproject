# Include this so that compile/run errors terminate script
#set -e

# Place all necessary fortran programs/modules into a single variable
myprogramfiles="command_line.f90 create_axis.f90 initialize.f90 potential.f90 motion.f90 write_netcdf.f90 main.f90"

# Name of the compiled file
outfile="electron"

# Name of compiler
fc=gfortran

# We use nf-config to grab the compile and link flags, backticks run command and grab output
fflags=`nf-config --fflags`
flibs=`nf-config --flibs`

# Command to compile program/modules along with necessary flags
$fc -g -std=f2008 -Wall $fflags $myprogramfiles $flibs -o $outfile

# Terminal command to run the simulation
./electron nx=100 ny=100 problem=single

# Plotting the visualizations of the electron path and potential
python3 netcdf_plotting.py
