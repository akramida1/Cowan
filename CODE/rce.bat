@echo off
cls
echo .
echo                              Program RCE Mod 20kd
echo                  Least-squares fitting of atomic energy levels
echo                                   ÜÜÜÜÜÜÜÜÜ
echo                                   Û START Û
echo                                   ßßßßßßßßß
rem rce --ground puts the calculated ground level of OUTGINE at zero
rem first (add_eav --ground).  A warning from it, as for INVECT not 0,
rem is no reason to stop; an error is, and RCE is then not run.
if /i "%~1" == "--ground" (
  call "%~dps0add_eav.bat" --ground
  if errorlevel 1 exit /b 1
)
call mark_time.bat %~dps0rce20kd.exe
if not exist RCEINP.HF copy RCEINP RCEINP.HF > nul
