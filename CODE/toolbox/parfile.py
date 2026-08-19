"""A parameter file in the author's own layout.

    <value>    ; <key>  <description of the parameter and its options>

The value is everything to the left of the first semicolon; the key is
the first word to the right of it, and the rest of the line is prose.
A line whose first non-blank character is a semicolon is a comment, and
so is a blank line, which is what lets a long description run on over
as many lines as it needs.  Order does not matter and a key that is
left out takes its default.

WHY NOT configparser, JSON OR TOML

The stock formats put the name first and the description nowhere: an
INI or TOML file carries `num_trials = 100` and the explanation has to
live in a comment that nothing checks, or in a manual.  This layout
puts the value where the eye finds it -- one column, all of them
aligned -- and makes the file its own documentation, which is the same
choice Cowan's fixed-form decks make and the one the author has been
using for his own codes.  The cost is this module, which is small, and
the gain is that `template()` can write a fresh, fully commented file
from the same table the reader validates against, so the two cannot
drift apart.

What this module adds over a hand-rolled reader is the checking: an
unknown key, a repeated key, a value that is not a number, and a number
outside its documented range are each named with the line they are on
and stop the run.  A parameter file is read once at the start of a
calculation that then takes hours; a typo in it must not be found by
the results being quietly wrong.
"""

from dataclasses import dataclass
from pathlib import Path

from . import files

#: The kinds a value may have.
INT, FLOAT, BOOL, TEXT = 'int', 'float', 'bool', 'text'

_TRUE = ('1', 'y', 'yes', 't', 'true', 'on')
_FALSE = ('0', 'n', 'no', 'f', 'false', 'off')

#: How wide the value column is in a file this module writes.
_COLUMN = 11


class ParFileError(Exception):
    """A parameter file that cannot be used as it stands."""


@dataclass(frozen=True)
class Spec:
    """One parameter: what it is called, what it may be, what it means."""

    key: str
    kind: str
    default: object
    description: str                # prose, wrapped by `template()`
    minimum: object = None
    maximum: object = None
    choices: tuple = ()

    def show(self, value):
        """The value as this module writes it back out."""
        if value is None:
            # A parameter with no default -- one that has to be filled
            # in before the program can run -- is written blank, and
            # reads back as None so the program can say what is missing.
            return ''
        if self.kind == BOOL:
            return '1' if value else '0'
        if self.kind == TEXT:
            return '' if value is None else str(value)
        if self.kind == INT:
            return '%d' % value
        return _trim(value)


def _trim(value):
    """A float in the shortest form that reads back as the same number."""
    text = repr(float(value))
    return text[:-2] if text.endswith('.0') else text


# --------------------------------------------------------------------------
# finding the file
# --------------------------------------------------------------------------

def find(name, folder='.', named=None, shipped=True):
    """The parameter file to use, or None when there is none.

    The order is the one every program in this package wants: a file
    named on the command line, then one in the working directory, then
    -- for a program that ships a filled-in default -- the one beside
    the code.  A `named` file that is not there is an error, because
    the user asked for that file and no other; nothing else is, so the
    caller decides what "none" means.  `conv_out` passes
    `shipped=False`: it must not ship a par file, or the file would
    never be absent and the sample would never be written.
    """
    if named:
        path = Path(named)
        if not path.is_file():
            raise ParFileError('%s is not there' % path)
        return path
    local = files.resolve(name, folder, required=False)
    if local is not None:
        return local
    if shipped:
        beside = files.code_dir() / name
        if beside.is_file():
            return beside
    return None


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def read(path, specs):
    """`{key: value}` for every spec, from the file at `path`.

    Every key of `specs` is present in the result whether or not the
    file mentions it; one the file leaves out takes its default.
    """
    table = {spec.key.lower(): spec for spec in specs}
    if len(table) != len(specs):
        raise ParFileError('two parameters of this program share a name')
    out = {}
    seen = {}
    for number, line in _lines(path):
        value, _, comment = line.partition(';')
        words = comment.split()
        key = words[0].lower() if words else ''
        if not key:
            raise ParFileError(
                '%s line %d: a value with no key.  The key is the first word '
                'after the semicolon, as in "100   ; num_trials  ..."'
                % (path, number))
        spec = table.get(key)
        if spec is None:
            raise ParFileError(
                '%s line %d: %r is not a parameter of this program.  It '
                'takes %s.' % (path, number, words[0],
                               ', '.join(sorted(table))))
        if key in seen:
            raise ParFileError('%s line %d: %s was already given on line %d'
                               % (path, number, spec.key, seen[key]))
        seen[key] = number
        out[spec.key] = _value(spec, value.strip(), path, number)
    for spec in specs:
        out.setdefault(spec.key, spec.default)
    return out


def _lines(path):
    """(line number, text) of every line that is not a comment."""
    with open(path, 'r', errors='replace') as fh:
        for number, raw in enumerate(fh, 1):
            line = raw.rstrip('\n').rstrip('\r')
            if not line.strip() or line.lstrip().startswith(';'):
                continue
            yield number, line


def _value(spec, text, path, number):
    """One field, converted and range-checked, or a complaint."""
    if not text:
        # A blank value is "use the default", which is how a parameter
        # is turned off without deleting the line that documents it.
        return spec.default
    if spec.kind == TEXT:
        if spec.choices and text.lower() not in [c.lower()
                                                 for c in spec.choices]:
            raise _bad(path, number, spec, text,
                       'one of ' + ', '.join(spec.choices))
        return text
    if spec.kind == BOOL:
        if text.lower() in _TRUE:
            return True
        if text.lower() in _FALSE:
            return False
        raise _bad(path, number, spec, text, '0 or 1')
    try:
        value = int(text) if spec.kind == INT else float(text)
    except ValueError:
        raise _bad(path, number, spec, text,
                   'a whole number' if spec.kind == INT else 'a number')
    if spec.minimum is not None and value < spec.minimum:
        raise _bad(path, number, spec, text, 'at least %s' % spec.minimum)
    if spec.maximum is not None and value > spec.maximum:
        raise _bad(path, number, spec, text, 'at most %s' % spec.maximum)
    return value


def _bad(path, number, spec, text, wanted):
    return ParFileError('%s line %d: %s = %r, but it must be %s'
                        % (path, number, spec.key, text, wanted))


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------

def template(specs, values=None, title=(), width=78):
    """The lines of a complete, commented parameter file.

    `values` supplies the values to write; anything it leaves out is
    written at its default.  The description of each parameter is
    wrapped into comment lines beneath it, so a file written by this
    function is a file this module can read back.
    """
    values = values or {}
    out = list(_banner(title, width))
    indent = ';' + ' ' * (_COLUMN + 2)
    for spec in specs:
        out.append('')
        value = values.get(spec.key, spec.default)
        out.append('%s; %s' % (spec.show(value).ljust(_COLUMN), spec.key))
        for piece in _wrap(spec.description, width - len(indent)):
            out.append((indent + piece).rstrip())
    return out


def _banner(title, width):
    yield ';' + '-' * (width - 1)
    for line in title:
        yield ('; ' + line).rstrip()
    if title:
        yield ';'
    for line in _wrap(
            'One parameter per line, as  <value>  ; <key>  <description>.  '
            'The value is everything left of the first semicolon and the key '
            'is the first word after it.  A blank value means "use the '
            'default".  A line starting with a semicolon is a comment, and '
            'so is a blank line.  Order does not matter.', width - 3):
        yield ('; ' + line).rstrip()
    yield ';' + '-' * (width - 1)


def _paragraphs(text):
    """Split prose into paragraphs on blank lines, ignoring the line breaks
    a triple-quoted string in the source happens to have in it."""
    out, current = [], []
    for line in text.split('\n'):
        if line.strip():
            current.append(line.strip())
        elif current:
            out.append(' '.join(current))
            current = []
    if current:
        out.append(' '.join(current))
    return out


def _wrap(text, width):
    """Wrap prose to `width`, keeping a blank line between paragraphs."""
    width = max(width, 24)
    out = []
    for paragraph in _paragraphs(text):
        if out:
            out.append('')
        line = ''
        for word in paragraph.split():
            if line and len(line) + 1 + len(word) > width:
                out.append(line)
                line = word
            else:
                line = '%s %s' % (line, word) if line else word
        if line:
            out.append(line)
    return out


def write(path, specs, values=None, title=(), ending='\n'):
    """Write a fresh parameter file; an existing one is never replaced."""
    path = Path(path)
    if path.exists():
        raise ParFileError('%s is already there; delete or rename it first'
                           % path)
    with open(path, 'w', newline='') as fh:
        fh.write(''.join(line + ending
                         for line in template(specs, values, title)))
    return path
