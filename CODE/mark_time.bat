@echo off
echo Started at                                   %time%, %date%
setlocal
set _time_ms=0
call :get_time
echo .
%1 %2 %3
echo .
call :get_time
set /A _hr=_time_ms / 360000
set /A _rest=_time_ms - _hr * 360000
set /A _min=_rest / 6000
set /A _rest=_rest - _min * 6000
set /A _sec=_rest / 100
set /A _ms= _rest - _sec * 100
echo Finished at                                  %time%, %date%
echo Time spent:                                  %_hr% hr %_min% min %_sec%.%_ms% sec
endlocal
goto :EOF

:get_time
for /f "tokens=1,2,3,4 delims=:." %%a in ("%time%") do (
  set /A _hr=%%a
  set /A _min=1%%b-100
  set /A _sec=1%%c-100
  set /A _ms=1%%d-100
  set /A _time_ms=_hr * 360000 + _min * 6000 + _sec * 100 + _ms - _time_ms
)
goto :EOF
