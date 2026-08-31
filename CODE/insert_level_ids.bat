@echo off
rem Put unique level identifiers into the tables conv_out writes.
setlocal
py "%~dps0insert_level_ids.py" %*
set _rc=%errorlevel%
if %_rc% == 9009 echo You do not have Python in your PATH.
rem A false `if` has just cleared ERRORLEVEL, which is why the exit code
rem was saved first and is handed back here rather than fallen out of.
endlocal & exit /b %_rc%
