@echo off
cls
echo .
echo .
echo .
echo ************************************************************************
echo *    This batch file creates the derived data files RCG reads: the     *
echo *      CFP binary decks FOR072, FOR073 and FOR074, and the LS-to-jj    *
echo *     deck LSJJ.DAT.  You need to run it once after installing the     *
echo *     package, and again if you re-compiled RCG with any compiler      *
echo *     different from Salford (Silverfrost) F77.  It uses the default   *
echo *                  CFP input deck (file ING11.CFP).                    *
echo ************************************************************************
rem
rem  Usage:  make_cfp.bat [target directory]
rem
rem  With no argument the decks are written beside this file, which is
rem  where RCG looks for them in an ordinary installation.  An argument
rem  writes them somewhere else instead, and SENIOR is copied there with
rem  them, so that the directory is a complete set of the permanent files
rem  RCG opens through COWAN.CFG.  That is for a run that must not write
rem  into the package directory -- the test suite is one, having no decks
rem  in a source tree and no licence to create any there.  RCG itself is
rem  always the one beside this file.
rem
rem  All four decks are derived files and none of them is distributed:
rem  the CFP decks and LSJJ.DAT alike are unformatted files whose record
rem  layout is the compiler's and whose numbers are the machine's, so
rem  they are built where they will be read.  LSJJ.DAT is built by
rem  lsjj\mkdeck.py and needs numpy; if numpy is absent this script says
rem  so and carries on, because the CFP decks are the primary job and
rem  everything except true jj coupling works without LSJJ.DAT.
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
:Sysfiles
if exist ING11 @echo Backing up existing ING11. to ING11.BAK
if exist ING11 copy ING11 ING11.BAK > nul
if exist ING11 del ING11 > nul
rem copy ING11.CFP ING11 > nul
echo                               CREATING ING11 . . .
for /F "skip=8 tokens=1* delims=;" %%a in (%SYSS%ING11.CFP) do @echo %%a>>ING11
rem  LSJJ.DAT is built from the ING11 just written and from SENIOR, and
rem  needs nothing else -- RCG is not involved, so it goes here, ahead of
rem  the RCG run, and a failure below only costs true jj coupling.
rem  cowan_cfp.py takes both inputs from these two variables when they
rem  are set, which is why no file has to be moved to make this work.
echo                             CREATING LSJJ.DAT . . .
set "LSJJ_CFP=%CD%\ING11"
set "LSJJ_SENIOR=%SYS%SENIOR"
set "LSJJOK="
rem  Three outcomes are distinguished, and none of them is fatal: built,
rem  no numpy, and anything else.  ERRORLEVEL is saved before the first
rem  `if` because a false `if` clears it.
py -c "import numpy" 2>nul
set "_rc=%errorlevel%"
if "%_rc%"=="9009" goto NoPython
if not "%_rc%"=="0" goto NoNumpy
py "%SYS%lsjj\mkdeck.py" "%OUT%LSJJ.DAT"
set "_rc=%errorlevel%"
if not "%_rc%"=="0" goto LsjjErr
set "LSJJOK=1"
goto Cfp
:NoPython
echo .
echo WARNING: you do not have Python in your PATH, so LSJJ.DAT was not
echo          rebuilt.  True jj coupling (column 38 of the ING11 main
echo          card set to 8) will not be available; everything else in
echo          the package works.  Install Python 3, then run make_cfp
echo          again.
echo .
goto LsjjFallback
:NoNumpy
echo .
echo WARNING: numpy is not installed, so LSJJ.DAT was not rebuilt.  True
echo          jj coupling (column 38 of the ING11 main card set to 8)
echo          will not be available; everything else in the package
echo          works.  Install it with
echo .
echo              pip install numpy
echo .
echo          and run make_cfp again.
echo .
goto LsjjFallback
:LsjjErr
echo .
echo WARNING: lsjj\mkdeck.py failed (exit code %_rc%; its own message is
echo          above), so LSJJ.DAT was not rebuilt.  True jj coupling
echo          will not be available; everything else in the package
echo          works.  The CFP decks are built below regardless.
echo .
:LsjjFallback
rem  The deck could not be built.  If we are writing somewhere other than
rem  the package directory and a deck is sitting in the package
rem  directory, that one is at least the right file for this machine, so
rem  fall back to copying it.
if /I "%OUT%"=="%SYS%" goto Cfp
if exist "%OUT%LSJJ.DAT" goto Cfp
if not exist "%SYS%LSJJ.DAT" goto Cfp
echo          Copying the existing "%SYS%LSJJ.DAT" instead.
copy "%SYS%LSJJ.DAT" "%OUT%" > nul
:Cfp
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
if defined LSJJOK echo       The LS-to-jj deck LSJJ.DAT was rebuilt there as well.
if not defined LSJJOK echo       LSJJ.DAT was NOT rebuilt -- see the warning above.
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
