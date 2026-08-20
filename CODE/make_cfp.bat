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
rem  Usage:  make_cfp.bat [target directory]
rem
rem  With no argument the decks are written beside this file, which is
rem  where RCG looks for them in an ordinary installation.  An argument
rem  writes them somewhere else instead, and SENIOR and LSJJ.DAT are
rem  copied there with them, so that the directory is a complete set of
rem  the permanent files RCG opens through COWAN.CFG.  That is for a run
rem  that must not write into the package directory -- the test suite is
rem  one, having no decks in a source tree and no licence to create any
rem  there.  RCG itself is always the one beside this file.
rem
setlocal
rem  SYS and OUT are the two directories, and each is kept twice: the
rem  name as given, for the messages, and the 8.3 name, for the two
rem  places that cannot carry a blank -- the unquoted for /F below, and
rem  the path RCG reads from COWAN.CFG and cuts at the first blank.
set "SYS=%~dp0"
set "SYSS=%~dps0"
set "OUT=%SYS%"
set "OUTS=%SYSS%"
if "%~1"=="" goto Ready
if not exist "%~f1" mkdir "%~f1"
if not exist "%~f1" goto Err3
set "OUT=%~f1\"
set "OUTS=%~fs1\"
:Ready
if not exist "%SYS%ING11.CFP" goto Err1
if /I "%OUT%"=="%SYS%" goto Sysfiles
if not exist "%OUT%SENIOR" copy "%SYS%SENIOR" "%OUT%" > nul
if not exist "%OUT%LSJJ.DAT" copy "%SYS%LSJJ.DAT" "%OUT%" > nul
:Sysfiles
if exist ING11 @echo Backing up existing ING11. to ING11.BAK
if exist ING11 copy ING11 ING11.BAK > nul
if exist ING11 del ING11 > nul
rem copy ING11.CFP ING11 > nul
echo                               CREATING ING11 . . .
for /F "skip=8 tokens=1* delims=;" %%a in (%SYSS%ING11.CFP) do @echo %%a>>ING11
if exist "%OUT%FOR072" del "%OUT%FOR07?" >nul
echo %OUTS% > COWAN.CFG
rem goto End
echo                                RUNNING RCG . . .
call mark_time.bat "%SYS%rcg11kd.exe"
rem del COWAN.CFG > nul
if exist FOR033 del FOR033 > nul
if exist TAPE2E del TAPE2E > nul
if exist OUTGINE del OUTGINE > nul
rem del ING11 > nul
if not exist "%OUT%FOR072" goto Err2
if not exist "%OUT%FOR073" goto Err2
if not exist "%OUT%FOR074" goto Err2
echo Done. New CFP binary files FOR072, FOR073, and FOR074 are created in "%OUT%" directory.
goto End
:Err1
echo Error: file "%SYS%ING11.CFP" is not found
goto End
:Err2
echo Files FOR072, FOR073 and FOR074 were not created in "%OUT%" due to execution error.
goto End
:Err3
echo Error: directory "%~f1" does not exist and could not be created
:End
endlocal
