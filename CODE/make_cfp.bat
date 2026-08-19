@echo off
cls
echo .
echo .
echo .
echo ************************************************************************
echo *      This batch file creates the new CFP binary files for RCG.       *
echo *     You need to run it if you re-compiled RCG with any compiler      *
echo *      different from Salford (Silverfrost) F77, to create new         *
echo *        FOR072, FOR073 and FOR074 files. It uses the default          *
echo *                  CFP input deck (file ING11.CFP).                    *
echo ************************************************************************
rem
setlocal
if not exist "%~dp0ing11.cfp" goto Err1
if exist ing11 @echo Backing up existing ING11. to ING11.BAK
if exist ing11 copy ing11 ing11.bak > nul
if exist ing11 del ing11 > nul
rem copy ing11.cfp ing11 > nul
echo                               CREATING ING11 . . .
for /F "skip=8 tokens=1* delims=;" %%a in (%~dps0ing11.cfp) do @echo %%a>>ing11
if exist "%~dp0for072" del "%~dp0for07?" >nul
echo %~dps0 > cowan.cfg
rem goto End
echo                                RUNNING RCG . . .
call mark_time.bat "%~dp0rcg11kd.exe"
rem del cowan.cfg > nul
if exist for033 del for033 > nul
if exist tape2e del tape2e > nul
if exist outgine del outgine > nul
rem del ing11 > nul
if not exist "%~dp0for072" goto Err2
if not exist "%~dp0for073" goto Err2
if not exist "%~dp0for074" goto Err2
echo Done. New CFP binary files for072, for073, and for074 are created in "%~dp0" directory.
goto End
:Err1
echo Error: file "%~dp0ing11.cfp" is not found
goto End
:Err2
echo Files FOR072, FOR073 and FOR074 were not created in "%~dp0" due to execution error.
:End
endlocal
