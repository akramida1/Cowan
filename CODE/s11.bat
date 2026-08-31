@echo off
rem Substitute fitted RCE parameters back into the ING11 deck.
rem The Perl this replaced is in perl_legacy\ (plan, phase 11).
setlocal
py "%~dps0s11.py" %*
set _rc=%errorlevel%
if %_rc% == 9009 echo You do not have Python in your PATH.
rem A false `if` has just cleared ERRORLEVEL, which is why the exit code
rem was saved first and is handed back here rather than fallen out of.
endlocal & exit /b %_rc%
