gfortran -g -O0 -fbacktrace -c param.f90 Functions.f90 cohini.f90 cohrun.f90
gfortran -g -O0 -fbacktrace boss.f90 param.o Functions.o cohini.o cohrun.o
a.exe 
