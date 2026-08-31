@echo off
rem  %1 is the source stem, one per call:  gfortran_O3_OMP rcg11kd
rem
rem  mkverres.bat writes %1_ver.o, a compiled VERSIONINFO resource, and is
rem  allowed to fail: it prints one line and produces nothing, and the link
rem  below then leaves the argument out.  Its header says why the resource
rem  is there and why CODE/unix/build.sh has no counterpart to it.
setlocal
call "%~dp0mkverres.bat" %1
set "VERRES="
if exist "%1_ver.o" set "VERRES=%1_ver.o"
gfortran.exe -std=legacy -fshort-enums -malign-double -ftracer -fno-backslash -fallow-argument-mismatch -fcray-pointer -O3 -fopenmp -fno-automatic -fallow-invalid-boz -o %1.exe %1.F %VERRES% 2>gfortran_out.txt
if exist "%1_ver.o" del "%1_ver.o"
endlocal
