@echo off
rem Write the table that translates level identifiers into IDEN2 numbers.
setlocal
py "%~dps0write_IDEN_level_ids.py" %*
set _rc=%errorlevel%
if %_rc% == 9009 echo You do not have Python in your PATH.
rem A false `if` has just cleared ERRORLEVEL, which is why the exit code
rem was saved first and is handed back here rather than fallen out of.
endlocal & exit /b %_rc%
