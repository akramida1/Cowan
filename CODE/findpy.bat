@echo off
rem Sets PY to the Python interpreter the utility wrappers run.  It is
rem called from inside their setlocal, so PY goes away with them, and it
rem has no setlocal of its own, which is what lets PY reach the caller.
rem
rem   1. A conda installation, found by `where conda`.  No conda variable
rem      outlives an activated session, so the test is for conda itself,
rem      and then:
rem        - in an activated session, that environment's python.exe
rem          (CONDA_PREFIX, which activation sets and deactivation clears);
rem        - otherwise the base installation's python.exe, which is one
rem          level above Scripts\conda.exe or condabin\conda.bat and two
rem          above Library\bin\conda.bat -- the three places a conda on
rem          the PATH is found;
rem        - otherwise plain `python`.
rem      The py launcher is passed over here even where it is installed:
rem      it runs the python.org interpreter it was installed with and
rem      knows nothing of conda, so its Python can lack numpy while the
rem      conda one has it.
rem   2. Otherwise the py launcher, which a python.org install provides.
rem   3. Otherwise plain `python`.
rem
rem A caller that finds exit code 9009 after "%PY%" ... has been told
rem that none of these exists, and says so.
set "PY="
where conda >nul 2>nul
if errorlevel 1 goto NoConda
if defined CONDA_PREFIX if exist "%CONDA_PREFIX%\python.exe" (
  set "PY=%CONDA_PREFIX%\python.exe"
  goto :eof
)
for /f "delims=" %%I in ('where conda 2^>nul') do (
  if not defined PY if exist "%%~dpI..\python.exe" set "PY=%%~dpI..\python.exe"
  if not defined PY if exist "%%~dpI..\..\python.exe" set "PY=%%~dpI..\..\python.exe"
)
rem  ..\ is resolved away, so that a message naming the interpreter
rem  names it plainly.
if defined PY for %%P in ("%PY%") do set "PY=%%~fP"
if not defined PY set "PY=python"
goto :eof
:NoConda
where py >nul 2>nul
if errorlevel 1 (set "PY=python") else (set "PY=py")
