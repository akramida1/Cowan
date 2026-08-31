@echo off
rem One coupling scheme per configuration, applied -> levels_prepub.txt.
rem There is no Perl ancestor; this tool is new (plan, phases 12-14).
setlocal
py "%~dps0make_levels_list_prepub.py" %*
set _rc=%errorlevel%
if %_rc% == 9009 echo You do not have Python in your PATH.
rem A false `if` has just cleared ERRORLEVEL, which is why the exit code
rem was saved first and is handed back here rather than fallen out of.
endlocal & exit /b %_rc%
