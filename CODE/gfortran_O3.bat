rem set LIBPATH="C:\MinGW\mingw64\lib"
gfortran.exe -std=legacy -fshort-enums -malign-double -ftracer -fno-backslash -fallow-argument-mismatch -fcray-pointer -O3 -o %1.exe %1.F 2>gfortran_out.txt
