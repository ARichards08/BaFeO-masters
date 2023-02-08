#!/bin/bash
#-------------------------------------------------------------------------------
#      Sample gnuplot template for plotting magnetization m(T) curve
#
#                (c) R F L Evans 2016. All rights reserved.
#
#-------------------------------------------------------------------------------

# Fit
f(x)= a + b*x + c*x**2
a=1;b=1;c=1
fit f(x) 'mag_moments.txt' u 1:2 via a, b, c
a=1;b=1;c=1
fit f(x) 'mag_moments.txt' u 1:3 via a, b, c
a=1;b=1;c=1
fit f(x) 'mag_moments.txt' u 1:4 via a, b, c
a=1;b=1;c=1
fit f(x) 'mag_moments.txt' u 1:5 via a, b, c
a=1;b=1;c=1
fit f(x) 'mag_moments.txt' u 1:6 via a, b, c
#--------------------------------------------
# Plot graph with fit
#p "timeperN.txt" u 1:2 ls 102 w l title "Runtime", \
# f(x) ls 103 w l title "Cubic Fit", \
# g(x) ls 104 w l title "Quadratic Fit"

