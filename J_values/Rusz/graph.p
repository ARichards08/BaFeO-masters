#!/bin/bash
#-------------------------------------------------------------------------------
#      Sample gnuplot template for plotting magnetization m(T) curve
#
#                (c) R F L Evans 2016. All rights reserved.
#
#-------------------------------------------------------------------------------

# Fit
f(x)= a*x**2 +b*x +c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 1:2 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 3:4 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 5:6 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 7:8 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 9:10 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 11:12 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 13:14 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 15:16 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 17:18 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 19:20 via a, b, c
a=1;b=1;c=1
fit f(x) 'J_points.txt' u 21:22 via a, b, c

#--------------------------------------------
# Plot graph with fit
#p "timeperN.txt" u 1:2 ls 102 w l title "Runtime", \
# f(x) ls 103 w l title "Cubic Fit", \
# g(x) ls 104 w l title "Quadratic Fit"

