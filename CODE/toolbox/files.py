"""Locating the files of a working directory, on any platform.

The decks in WORK/ spell the same file both ways - ING11 and ing11 -
because Windows does not distinguish them.  Linux and macOS do.  Every
tool in this package therefore opens deck files through resolve() rather
than by name, so one deck runs unchanged on all three systems.

Fortran writes fixed-form text with trailing blanks and, on Windows, CRLF
line endings; read_lines() strips both so that column positions counted
from the FORMAT statements land where they should.
"""

from pathlib import Path

# The files of a run, in the spelling the Fortran itself uses.
IN36 = 'IN36'
IN2 = 'IN2'
ING11 = 'ING11'
OUTG11 = 'OUTG11'
OUT36 = 'OUT36'
TAPE2N = 'TAPE2N'
GENEALOGY = 'genealogy.txt'
JJKEYS = 'true-jj_keys.txt'


def resolve(name, folder='.', required=True):
    """Return the path of `name` in `folder`, ignoring case.

    An exact match wins; otherwise the directory is scanned once.  Two
    files differing only in case are an error rather than a coin toss.
    """
    folder = Path(folder)
    exact = folder / name
    if exact.is_file():
        return exact
    lowered = name.lower()
    hits = [p for p in folder.iterdir() if p.is_file() and p.name.lower() == lowered] \
        if folder.is_dir() else []
    if len(hits) > 1:
        names = ', '.join(sorted(p.name for p in hits))
        raise FileExistsError(
            'several files in %s differ only in case: %s' % (folder, names))
    if hits:
        return hits[0]
    if required:
        raise FileNotFoundError('%s not found in %s' % (name, folder))
    return None


def read_lines(path, width=None, strip=True):
    """Read a fixed-form text file as a list of lines, CR and blanks off.

    With `width` the lines are blank-padded to that many columns, so a
    reader may index any documented column without a length test.

    `strip=False` keeps the trailing blanks.  A reader wants them gone;
    a **rewriter** wants them, because a card it does not touch has to
    come back out as the bytes it came in as, and Fortran pads a record
    to its FORMAT's full width whether or not anything is in it.
    """
    with open(path, 'r', errors='replace') as fh:
        lines = fh.read().split('\n')
    if lines and lines[-1] == '':
        lines.pop()
    out = [ln.rstrip('\r') for ln in lines]
    if strip:
        out = [ln.rstrip() for ln in out]
    if width:
        out = [ln.ljust(width) for ln in out]
    return out


def line_ending(path, folder=None):
    """The line terminator a file uses, so that a rewrite keeps it.

    RCG writes CRLF on Windows and LF elsewhere, and Fortran reads
    either; but a copy that changes it differs from its original on
    every line, which is exactly what a rewritten deck must not do.
    """
    if folder is not None:
        path = resolve(path, folder)
    with open(path, 'rb') as fh:
        head = fh.read(4096)
    return '\r\n' if b'\r\n' in head else '\n'


def write_deck(path, lines, ending='\n', backup=None):
    """Write `lines` to `path`, keeping a copy of what was there.

    The Perl utilities shell out to `del`, `ren` and `copy` for this,
    which is one of the two things that tie them to Windows.  The new
    file is written beside the old one and moved into place, so a run
    that fails partway leaves the deck it was given intact.

    `backup` is the name the old file is kept under.  An existing backup
    of that name is replaced, which is what the Perl does; the deck
    itself is never removed before the new one is complete.
    """
    path = Path(path)
    temporary = path.with_name(path.name + '.new')
    with open(temporary, 'w', newline='') as fh:
        fh.write(''.join(line + ending for line in lines))
    if backup is not None and path.exists():
        kept = path.with_name(backup)
        if kept.exists():
            kept.unlink()
        path.replace(kept)
    temporary.replace(path)
    return path


def code_dir():
    """The directory this package is installed in (CODE/, or its copy)."""
    return Path(__file__).resolve().parent.parent
