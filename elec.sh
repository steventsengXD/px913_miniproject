# Include this so that compile/run errors terminate script
#set -e

# Place all necessary fortran programs/modules into a single variable
myprogramfiles="command_line.f90 create_axis.f90 helpers.f90 write_netcdf.f90 main.f90"

# Name of the compiled file
outfile="electron"

# Name of compiler
fc=gfortran

# We use nf-config to grab the compile and link flags, backticks run command and grab output
fflags=`nf-config --fflags`
flibs=`nf-config --flibs`

# Command to compile program/modules along with necessary flags
$fc -g -std=f2008 -Wall $fflags $myprogramfiles $flibs -o $outfile
