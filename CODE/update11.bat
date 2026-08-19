copy outgine outgine.bak
call s11.bat %1
call rcg.bat
copy outgine.bak OUTGINE
del outgine.bak
call printout.bat levels.txt
