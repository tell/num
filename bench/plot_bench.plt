set terminal postscript eps enhanced monochrome blacktext \
   dashed dashlength 2.0 linewidth 4.0 defaultplex \
   palfuncparam 2000,0.003 \
   butt "Helvetica" 14

datname="kronecker-jacobi.bench.dat"

##

set output "NTZ.eps"
set xlabel "bitlength [bit]"
set ylabel "clock cycles [clk]"
plot datname ind 0 using 1:2 title "mpz\\_scan1" with lines, \
     datname ind 0 using 1:3 title "NTZ" with lines

set output "NTZ_ratio.eps"
set ylabel "ratio"
plot datname ind 0 using 1:($3/$2) title "NTZ / mpz\\_scan1" with lines

##

set output "shr.eps"
set ylabel "clock cycles [clk]"
plot datname ind 1 using 1:2 title "mpz\\_tdiv\\_q\\_2exp" with lines, \
     datname ind 1 using 1:3 title "shr" with lines, \
     datname ind 5 using 1:3 title "shr opt 1" with lines, \
     datname ind 9 using 1:3 title "shr opt 4" with lines

##

set output "sub.eps"
set ylabel "clock cycles [clk]"
plot datname ind 2 using 1:2 title "mpz\\_sub" with lines, \
     datname ind 2 using 1:3 title "sub" with lines, \
     datname ind 6 using 1:3 title "sub opt 1" with lines, \
     datname ind 10 using 1:3 title "sub opt 4" with lines

##

set output "kronecker-binary.eps"
set ylabel "clock cycles [clk]"
plot datname ind 3 using 1:2 title "mpz\\_kronecker" with lines, \
     datname ind 3 using 1:3 title "kronecker" with lines, \
     datname ind 7 using 1:3 title "kronecker using opt" with lines, \
     datname ind 11 using 1:3 title "kronecker using opt 4" with lines

##

# not yet
# set output "NTZ_opt.eps"

####

set xlabel "bitlength [bit]"
set ylabel "ratio"

set output "shr_ratio.eps"
plot datname ind 1 using 1:($3/$2) title "shr / mpz\\_tdiv\\_q\\_2exp" with lines, \
     datname ind 5 using 1:($3/$2) title "shr opt 1 / mpz\\_tdiv\\_q\\_2exp" with lines, \
     datname ind 9 using 1:($3/$2) title "shr opt 4 / mpz\\_tdiv\\_q\\_2exp" with lines

set output "sub_ratio.eps"
plot datname ind 2 using 1:($3/$2) title "sub / mpz\\_sub" with lines, \
     datname ind 6 using 1:($3/$2) title "sub opt 1 / mpz\\_sub" with lines, \
     datname ind 10 using 1:($3/$2) title "sub opt 4 / mpz\\_sub" with lines

set output "kronecker-binary_ratio.eps"
plot datname ind 3 using 1:($3/$2) title "kronecker / mpz\\_kronecker" with lines, \
     datname ind 7 using 1:($3/$2) title "kronecker using opt / mpz\\_kronecker" with lines, \
     datname ind 11 using 1:($3/$2) title "kronecker using opt 4 / mpz\\_kronecker" with lines
