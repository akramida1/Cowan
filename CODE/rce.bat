@echo off
cls
echo .
echo                              Program RCE Mod 20kd
echo                  Least-squares fitting of atomic energy levels
echo                                   ÜÜÜÜÜÜÜÜÜ
echo                                   Û START Û
echo                                   ßßßßßßßßß
rem Every file RCE writes must be writable before it starts.  gfortran
rem opens a file it cannot write -- one held open by Excel, or marked
rem read-only -- for reading instead, without a word, and fails only at
rem the first WRITE.  For LEVELS2 and LEVELS3 that is the very end of the
rem run, after LEVELS1, OUTE and RCEOUT are complete: the run looks
rem finished and leaves the two sorted lists of the previous run behind.
rem The test appends nothing, so it changes no file, and creates none.
rem RCEINP is not tested: RCE writes it only when INVECT is 0.
rem (`call set` re-expands _locked as each name is added to it.)
set "_locked="
for %%F in (OUTE PARVALS RCEOUT LEVELS1 LEVELS2 LEVELS3 LEVELS1J LEVELS2J LEVELS3J) do (
  if exist %%F ((call ) >>%%F) 2>nul || call set "_locked=%%_locked%% %%F"
)
if defined _locked (
  echo RCE not run: cannot write%_locked%.
  echo Close the program that has the file open, or clear its read-only
  echo attribute, and run rce again.
  set "_locked="
  exit /b 1
)
rem rce --ground puts the calculated ground level of OUTGINE at zero
rem first (add_eav --ground).  A warning from it, as for INVECT not 0,
rem is no reason to stop; an error is, and RCE is then not run.
if /i "%~1" == "--ground" (
  call "%~dps0add_eav.bat" --ground
  if errorlevel 1 exit /b 1
)
call mark_time.bat %~dps0rce20kd.exe
if not exist RCEINP.HF copy RCEINP RCEINP.HF > nul
