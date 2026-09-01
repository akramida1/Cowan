@echo off
rem Carry a level optimization's result back into the IDEN2 working files.
setlocal
py "%~dps0update_IDEN.py" %*
set _rc=%errorlevel%
if %_rc% == 9009 echo You do not have Python in your PATH.
rem A false `if` has just cleared ERRORLEVEL, which is why the exit code
rem was saved first and is handed back here rather than fallen out of.
endlocal & exit /b %_rc%
