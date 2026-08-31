@echo off
rem Level and line lists from an RCG/RCE run, with an intensity model.
rem Settings may also be kept in a conv_out.par in the working
rem directory; run it with no arguments to be given a commented one.
rem The Perl this replaced is in perl_legacy\ (plan, phase 11).
setlocal
py "%~dps0conv_out.py" %*
set _rc=%errorlevel%
if %_rc% == 9009 echo You do not have Python in your PATH.
rem A false `if` has just cleared ERRORLEVEL, which is why the exit code
rem was saved first and is handed back here rather than fallen out of.
endlocal & exit /b %_rc%
