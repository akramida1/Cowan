"""The LS->jj coupling algebra, and the generator of the LSJJ.DAT deck.

The four modules here are the production half of the work whose record is
in FOR/LSJJ/README.txt.  They live under CODE/ rather than beside that
record because make_cfp rebuilds LSJJ.DAT on the installed machine, and an
installed tree is CODE/ alone -- FOR/ is not copied by install.sh, nor by
the author's habit of copying CODE/ to the directory on the Windows PATH.

    lsjj.py       the LS and jj shell algebra
    cowan_cfp.py  Cowan's own LS cfp tables, read from the ING11 card deck
    recur.py      the jj cfp recursion
    mkdeck.py     writes LSJJ.DAT; run as a script, by make_cfp

They import one another as flat top-level modules (`import lsjj`), which is
what running mkdeck.py by path gives them: the interpreter puts this
directory first on sys.path.  This package exists so that the directory is
installed and found, not to be imported as one.  The validation scripts in
FOR/LSJJ/ reach these modules by adding this directory to sys.path.
"""
