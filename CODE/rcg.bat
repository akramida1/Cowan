@echo off
cls
echo %~dps0 > COWAN.CFG
echo .
echo                             Program RCG Mod 11kd
echo       computes angular matrix elements, energy levels, and atomic spectra
echo                                  ÜÜÜÜÜÜÜÜÜ
echo                                  Û START Û
echo                                  ßßßßßßßßß
call mark_time.bat "%~dp0rcg11kd.exe"
del COWAN.CFG > nul
