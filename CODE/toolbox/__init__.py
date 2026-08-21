"""Shared readers and formatters for R. D. Cowan's atomic-structure package.

The modules here replace the parsing that used to live, three times over,
inside the Perl utilities in CODE/.  Nothing in this package prints a
result by itself; the scripts one level up do that.

The Fortran FORMAT statements are the specification for every file read
here.  Each reader names the FOR/*.f line of the FORMAT it mirrors, so a
change in the Fortran can be traced to the code that has to follow it.
"""

__all__ = ['ai', 'asd', 'case', 'cmpnum', 'designate', 'files', 'fortfmt', 'fundconst',
           'genealogy', 'iden', 'in36', 'ing11', 'jjkeys', 'levels', 'outg11',
           'parnames', 'radiate', 'rce', 'spectrum', 'vacair']
