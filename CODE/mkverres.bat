@echo off
rem  Generate and compile the Windows VERSIONINFO resource for one of the
rem  four programs.
rem
rem      mkverres RCN2K      ->  %TEMP%\RCN2K_ver.rc  ->  .\RCN2K_ver.o
rem
rem  It is called by gfortran_O3.bat and gfortran_O3_OMP.bat, which link
rem  the .o beside the Fortran object and then delete it.  gfortran_dbg.bat
rem  does not call it: a debug build is never distributed and the resource
rem  would only be noise there.
rem
rem  WHY THE RESOURCE EXISTS.  Norton's Win64:Evo-gen heuristic deletes a
rem  freshly built RCN2K.exe every time -- measured, six trials out of six,
rem  and on binaries it had never seen, so it is a verdict on file content
rem  and not on prevalence or on the missing signature.  Two controls
rem  separate the cause cleanly: a resource carrying no version information
rem  at all survives, while an equal-sized padding of ordinary data with no
rem  .rsrc section is deleted.  The discriminating feature is the presence
rem  of a resource directory.  A PE without one is unusual for legitimate
rem  Windows software and common in malware, so the honest fix is to put
rem  there what the section is for.  All four programs get it, not only
rem  RCN2K: it is the same one line in each build, the information is true
rem  of each program, and it removes the feature before some future edit
rem  pushes another binary over the same threshold.
rem
rem  NOTHING HERE IS ALLOWED TO BREAK A BUILD.  No windres, no version.txt,
rem  an unknown stem -- each prints one line and exits 0 with no .o, and
rem  the caller links without it.
rem
rem  UNIX HAS NO COUNTERPART.  ELF carries no VERSIONINFO, so
rem  CODE/unix/build.sh is untouched by any of this.  It is one of the few
rem  places where the "Windows and Unix move together" rule of CLAUDE.md
rem  section 8 does not apply, because there is no other half to move.
rem
setlocal
if "%~1"=="" goto NoStem
set "SYS=%~dp0"
if not exist "%SYS%version.txt" goto NoFile
where windres >nul 2>nul
if errorlevel 1 goto NoRes
rem  version.txt is key=value, one per line, # for a comment.  tokens=1,*
rem  keeps the value whole, spaces, commas and full stops included.
set "V_VERSION="
set "V_COMMA_VERSION="
set "V_COMPANY="
set "V_PRODUCT="
set "V_COPYRIGHT="
for /F "usebackq eol=# tokens=1,* delims==" %%a in ("%SYS%version.txt") do set "V_%%a=%%b"
if not defined V_VERSION goto NoFile
if not defined V_COMMA_VERSION goto NoFile
rem  Four programs, and there always will be four, so a four-arm if is
rem  the whole mapping.  The stem is matched case-insensitively because
rem  the shipped names are not all capitalized the same way.
set "DESC="
if /I "%~1"=="RCN36K"  set "DESC=RCN -- Hartree-Fock atomic radial wavefunctions"
if /I "%~1"=="RCN2K"   set "DESC=RCN2 -- configuration-interaction integrals"
if /I "%~1"=="rcg11kd" set "DESC=RCG -- atomic energy levels and transition arrays"
if /I "%~1"=="RCE20KD" set "DESC=RCE -- least-squares fitting of energy parameters"
if not defined DESC goto NoStem
set "RC=%TEMP%\%~1_ver.rc"
rem  The .rc is written outside every parenthesised block on purpose:
rem  echo and the .rc grammar disagree about ( ) < > & | ^ and %, and a
rem  line ending in ) inside an if(...) block is the classic failure.
rem  Note also that a space left before >> becomes part of the line, so
rem  every redirection below is written hard against the text.
> "%RC%" echo 1 VERSIONINFO
>>"%RC%" echo FILEVERSION %V_COMMA_VERSION%
>>"%RC%" echo PRODUCTVERSION %V_COMMA_VERSION%
>>"%RC%" echo FILEFLAGSMASK 0x3fL
>>"%RC%" echo FILEFLAGS 0x0L
>>"%RC%" echo FILEOS 0x40004L
>>"%RC%" echo FILETYPE 0x1L
>>"%RC%" echo FILESUBTYPE 0x0L
>>"%RC%" echo BEGIN
>>"%RC%" echo   BLOCK "StringFileInfo"
>>"%RC%" echo   BEGIN
>>"%RC%" echo     BLOCK "040904b0"
>>"%RC%" echo     BEGIN
>>"%RC%" echo       VALUE "CompanyName", "%V_COMPANY%"
>>"%RC%" echo       VALUE "FileDescription", "%DESC%"
>>"%RC%" echo       VALUE "FileVersion", "%V_VERSION%"
>>"%RC%" echo       VALUE "InternalName", "%~1"
>>"%RC%" echo       VALUE "LegalCopyright", "%V_COPYRIGHT%"
>>"%RC%" echo       VALUE "OriginalFilename", "%~1.exe"
>>"%RC%" echo       VALUE "ProductName", "%V_PRODUCT%"
>>"%RC%" echo       VALUE "ProductVersion", "%V_VERSION%"
>>"%RC%" echo     END
>>"%RC%" echo   END
>>"%RC%" echo   BLOCK "VarFileInfo"
>>"%RC%" echo   BEGIN
>>"%RC%" echo     VALUE "Translation", 0x409, 1200
>>"%RC%" echo   END
>>"%RC%" echo END
windres -i "%RC%" -O coff -o "%~1_ver.o"
if errorlevel 1 goto Failed
if not exist "%~1_ver.o" goto Failed
del "%RC%" >nul 2>nul
endlocal & exit /b 0
:NoStem
echo mkverres: no version resource for "%~1" -- building without one.
endlocal & exit /b 0
:NoFile
echo mkverres: version.txt is missing or incomplete -- building without a
echo           version resource.
endlocal & exit /b 0
:NoRes
echo mkverres: windres is not on the PATH -- building without a version
echo           resource.  It ships with the same MinGW as gfortran.
endlocal & exit /b 0
:Failed
echo mkverres: windres could not compile "%RC%" -- building without a
echo           version resource.
endlocal & exit /b 0
