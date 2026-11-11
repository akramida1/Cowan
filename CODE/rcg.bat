@echo off
cls
rem  RCG opens the CFP decks FOR072, FOR073 and FOR074 beside this file
rem  with STATUS=UNKNOWN, so a deck that make_cfp never built is created
rem  empty and the first read stops with "End of file" on FOR072 -- which
rem  reads as a corrupt deck, not a missing one.  Check before that file
rem  is made.  The size modifier on a for variable is blank for a file
rem  that does not exist and 0 for an empty one.
for %%F in ("%~dp0FOR072" "%~dp0FOR073" "%~dp0FOR074") do (
  if "%%~zF"=="" goto NoDecks
  if "%%~zF"=="0" goto NoDecks
)
echo %~dps0 > COWAN.CFG
echo .
echo                             Program RCG Mod 11kd
echo       computes angular matrix elements, energy levels, and atomic spectra
echo                                  ÜÜÜÜÜÜÜÜÜ
echo                                  Û START Û
echo                                  ßßßßßßßßß
call mark_time.bat "%~dp0rcg11kd.exe"
del COWAN.CFG > nul
goto :EOF
:NoDecks
echo .
echo ERROR: the CFP decks FOR072, FOR073 and FOR074 are missing or empty in
echo        "%~dp0"
echo        RCG cannot run without them.  Build them by running make_cfp once,
echo        in a directory with nothing in it you need: it replaces ING11 there
echo        and deletes TAPE2E and OUTGINE.
echo .
exit /b 1
