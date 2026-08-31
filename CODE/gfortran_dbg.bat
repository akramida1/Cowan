@echo off
rem Debug build: gfortran_dbg <source stem>, e.g. gfortran_dbg RCE20KD
rem
rem The flags are gfortran_O3_OMP.bat's, with -O0 -ggdb3 in place of -O3
rem and without -fopenmp: a debugger is far easier to follow serially, and
rem the codes are bit-identical threaded or not (CLAUDE.md 1).
rem
rem -std=legacy and -fallow-invalid-boz are what let Cowan's sources
rem compile at all, and -fno-automatic is kept even without -fopenmp
rem because it puts locals in static storage -- drop it and the numbers
rem change, so a debug build would no longer be debugging the same code.
rem
rem -fcheck=bounds is deliberately absent: EQUIVALENCE and the dummy
rem arrays this package is built on trip it on correct code.
gfortran.exe -std=legacy -fshort-enums -malign-double -ftracer -fno-backslash -fallow-argument-mismatch -fcray-pointer -fallow-invalid-boz -fno-automatic -O0 -ggdb3 -fbacktrace -o %1.exe %1.F 2>gfortran_out.txt
