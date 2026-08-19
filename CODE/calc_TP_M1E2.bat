@echo off
rem Transition probabilities for forbidden (M1, E2) lines.
rem The physics settings are in calc_TP_M1E2.par, not here.
rem The Perl this replaced is in perl_legacy\ (plan, phase 11).
setlocal
py "%~dps0calc_TP_M1E2.py" %*
set _rc=%errorlevel%
if %_rc% == 9009 echo You do not have Python in your PATH.
rem A false `if` has just cleared ERRORLEVEL, which is why the exit code
rem was saved first and is handed back here rather than fallen out of.
endlocal & exit /b %_rc%
