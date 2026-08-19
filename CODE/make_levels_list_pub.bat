@echo off
rem Designations unique within each (parity, J) -> levels_pub.txt.
rem There is no Perl ancestor; this tool is new (plan, phases 12-14).
setlocal
py "%~dps0make_levels_list_pub.py" %*
set _rc=%errorlevel%
if %_rc% == 9009 echo You do not have Python in your PATH.
rem A false `if` has just cleared ERRORLEVEL, which is why the exit code
rem was saved first and is handed back here rather than fallen out of.
endlocal & exit /b %_rc%
