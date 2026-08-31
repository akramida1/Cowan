copy OUTGINE OUTGINE.BAK
rem With nothing typed, s11 takes the number of parity sets from IN36
rem ("auto"): a deck of both parities gets both sets substituted without
rem the trailing 2 the Perl needed.  An explicit 2, or a folder, still
rem goes through.
if "%~1" == "" (call s11.bat auto) else (call s11.bat %*)
call rcg.bat
copy OUTGINE.BAK OUTGINE
del OUTGINE.BAK
call printout.bat levels.txt
