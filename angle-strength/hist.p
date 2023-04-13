#!/bin/bash
#-------------------------------------------------------------------------------
#      Sample gnuplot template for plotting magnetization m(T) curve
#
#                (c) R F L Evans 2016. All rights reserved.
#
#-------------------------------------------------------------------------------
#Set terminal
set terminal postscript portrait enhanced color dashed lw 1 "DejaVuSans" 12

# Load styles
load "nature.journal"
load "bluegold-sl.colour"

# Set output file base name
ofile="angle-hist"

# Set axis labels
set xlabel "Bond angle (degrees)"
# For y-label use italics for variables and offset from axis
set ylabel "Frequency Density" offset 0.5,0

# Set axis number formatting
set format y "%3.2f"

# Set '0' line as guide to the eye
set xzeroaxis

# Set key properties
set key top right maxrows 5 spacing 1.2

# Set eps file name by concatenation (. operator)
epsfile=ofile.".eps"
set output epsfile

set xrange [0:180]

binwidth=10
bin(x, width)=width*floor(x/width)

#--------------------------------------------
# Plot graph with fit
p \
"angle-strength.dat" u (bin($3, binwidth)):(1.0/binwidth) smooth freq with boxes notitle
