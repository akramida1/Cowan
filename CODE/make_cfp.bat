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
if not exist "%~dp0ING11.CFP" goto Err1
if exist ING11 @echo Backing up existing ING11. to ING11.BAK
if exist ING11 copy ING11 ING11.BAK > nul
if exist ING11 del ING11 > nul
rem copy ING11.CFP ING11 > nul
echo                               CREATING ING11 . . .
for /F "skip=8 tokens=1* delims=;" %%a in (%~dps0ING11.CFP) do @echo %%a>>ING11
if exist "%~dp0FOR072" del "%~dp0FOR07?" >nul
echo %~dps0 > COWAN.CFG
rem goto End
echo                                RUNNING RCG . . .
call mark_time.bat "%~dp0rcg11kd.exe"
rem del COWAN.CFG > nul
if exist FOR033 del FOR033 > nul
if exist TAPE2E del TAPE2E > nul
if exist OUTGINE del OUTGINE > nul
rem del ING11 > nul
if not exist "%~dp0FOR072" goto Err2
if not exist "%~dp0FOR073" goto Err2
if not exist "%~dp0FOR074" goto Err2
echo Done. New CFP binary files FOR072, FOR073, and FOR074 are created in "%~dp0" directory.
goto End
:Err1
echo Error: file "%~dp0ING11.CFP" is not found
goto End
:Err2
echo Files FOR072, FOR073 and FOR074 were not created in "%~dp0" due to execution error.
:End
endlocal
