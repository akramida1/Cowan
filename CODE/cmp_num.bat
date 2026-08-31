@echo off
rem Compare two Fortran outputs, or two run folders, within a tolerance.
rem The oracle for changes that are not meant to be bit-exact; `diff`
rem stays the test for the ones that are (CLAUDE.md 1).
setlocal
py "%~dps0cmp_num.py" %*
set _rc=%errorlevel%
if %_rc% == 9009 echo You do not have Python in your PATH.
rem A false `if` has just cleared ERRORLEVEL, which is why the exit code
rem was saved first and is handed back here rather than fallen out of.
endlocal & exit /b %_rc%
