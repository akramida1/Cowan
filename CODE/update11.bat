copy OUTGINE OUTGINE.BAK
call s11.bat %1
call rcg.bat
copy OUTGINE.BAK OUTGINE
del OUTGINE.BAK
call printout.bat levels.txt
