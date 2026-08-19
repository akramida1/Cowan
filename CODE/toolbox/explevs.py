"""Carrying observed energies from one RCE deck into another.

A least-squares fit needs every level it fits to be tied to the right
eigenvector, and the tie has to survive the calculation changing --
another configuration added, a parameter freed, a different scaling --
which moves the fitted energies and reorders the levels.  What does not
move is the **eigenvector**: a level is recognised by the configurations
and terms it is built from and by how much of each.

That is what this module does.  For each observed level of one deck it
finds the level of another whose composition is nearest, and reports the
ones it is not sure about rather than assigning them quietly.

The distance is `subst_exp_levs.pl:294-336`: over the components of both
vectors, the sum of squared differences of the percentages, with a
component present in one and not the other counting its own square.  It
is computed twice, once with the second vector's signs as they are and
once with them reversed, and the smaller is taken -- an eigenvector is
only defined up to its overall sign, and RCE does not fix one.  Two
vectors sharing a component above 70 in both are the same level and the
search stops there.

The percentages are what RCE prints: a signed amplitude percentage,
whose square over 100 is the intensity percentage a table shows.
"""

import re
from dataclasses import dataclass, field

from . import files, rce

# A component below this is not worth carrying: `subst_exp_levs.pl:110`
# keeps only those at or past 10, and stops at the first one that is not.
FLOOR = 10

# Both vectors holding this much of one basis state settle it
# (`subst_exp_levs.pl:308`).  70 in amplitude is 49% in intensity, and
# two components cannot both be that large.
CERTAIN = 70

# RCE prints five components on a level line and `subst_exp_levs.pl:74-86`
# reads all five.
COMPONENTS = 5

_L = 'SPDFGHIKLMNOPQRSTUV'


class ExpLevelsError(Exception):
    """Two decks that cannot be matched level for level."""


@dataclass
class Level:
    """One level of an RCE deck, as the matcher needs it."""

    j: str                          # as written: '3', '3/2'
    ee: float                       # observed energy, or the computed one
    ec: float                       # fitted energy
    observed: bool
    vector: dict = field(default_factory=dict)   # {'conf,term': percent}
    kset: int = 1
    position: int = 0               # its place within its J, 1-based
    line: int = -1                  # the line it was read from
    assigned: float = None          # the energy this level has been given
    match: object = None            # the level of the other deck it is tied to

    @property
    def taken(self):
        """Whether this level already carries an observed energy.

        A level the file itself marks as observed is as unavailable as
        one this run has just assigned: both already stand for an
        observation, and neither may take a second.  The Perl keeps the
        two in one flag by clearing its observed mark on assignment
        (`subst_exp_levs.pl:342`).
        """
        return self.observed or self.assigned is not None

    @property
    def energy(self):
        """The observed energy this level carries, assigned or its own."""
        return self.ee if self.assigned is None else self.assigned

    @property
    def leading(self):
        """(percent, 'conf,term') of the largest component."""
        if not self.vector:
            return (0, '')
        key = max(self.vector, key=lambda k: abs(self.vector[k]))
        return (self.vector[key], key)

    def sorted_vector(self):
        """The components, largest first."""
        return sorted(self.vector.items(), key=lambda kv: -abs(kv[1]))


class Deck:
    """One RCE file's levels, grouped by parity set and J."""

    def __init__(self, path):
        self.path = path
        self.sets = {}              # {kset: {J: [Level]}}

    def add(self, level):
        self.sets.setdefault(level.kset, {}).setdefault(level.j, []) \
            .append(level)
        level.position = len(self.sets[level.kset][level.j])
        return level

    def levels(self):
        for kset in sorted(self.sets):
            for j in sorted(self.sets[kset], key=_j_value):
                for level in self.sets[kset][j]:
                    yield level

    def group(self, kset, j):
        return self.sets.get(kset, {}).get(j, [])


def _j_value(text):
    """'3/2' -> 1.5, for sorting a J that is written as a fraction."""
    if '/' in text:
        top, _, bottom = text.partition('/')
        return float(top) / float(bottom)
    return float(text)


def _j_text(j):
    """A J as RCE writes it on a level line: '3', '3/2'."""
    doubled = int(round(2 * j))
    return '%d' % (doubled // 2) if doubled % 2 == 0 else '%d/2' % doubled


def _key(component):
    """The '<configuration>,<term>' a component is matched on.

    Runs of blanks inside the twelve-column term label are squeezed to
    one so that the key does not depend on how wide RCE happened to pad
    it -- the `"$9 $10"` of `subst_exp_levs.pl:74`.
    """
    return '%s,%s' % (component.configuration.strip(),
                      ' '.join(component.label.split()))


def _vector(components):
    """{key: percent} for the components worth keeping."""
    out = {}
    for i, component in enumerate(components[:COMPONENTS]):
        if i and abs(component.percent) < FLOOR:
            break
        out[_key(component)] = component.percent
    return out


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def read_deck(path):
    """Every level of an RCEINP or RCEOUT, tied to the line it is on."""
    deck = Deck(path)
    csets = rce.read_rceout(path)
    lines = files.read_lines(path, width=140)
    places = {}
    for index, kset, position in rce.level_lines(lines):
        places[(kset, position)] = index
    for cset in csets:
        for level in cset.levels:
            deck.add(Level(j=_j_text(level.j), ee=level.e_obs,
                           ec=level.e_calc, observed=level.observed,
                           vector=_vector(level.components), kset=cset.index,
                           line=places.get((cset.index, level.index), -1)))
    return deck


# A LOPT level list: ten tab-separated fields, of which the matcher
# wants the designation, J, the energy and the second uncertainty, and
# then the configuration and term of the last two columns
# (`subst_exp_levs.pl:181-183`).
_LOPT = re.compile(r'^([^\t]*\t){9}[^\t]*$')
KK = 0.001                          # LOPT writes cm-1, RCE 1000 cm-1


def read_lopt(path, most=None):
    """Every level of a LOPT output, in the shape the matcher wants.

    LOPT is the level-optimization program the observed energies come
    out of, so this is how a freshly optimized level list is carried
    into a fit.  A level whose second uncertainty is above `most` is
    left out, which is the `max_unc` of `subst_exp_levs.pl:187`.

    Parity comes from the term's own asterisk, and which parity set it
    is filed under is settled by the first level read, exactly as
    `:197-199` settles it -- a LOPT list does not say.

    A row with no J is a **term** rather than a level, and one level is
    made for each J the term has, which is `:206-221`.

    There is no LOPT list among the reference cases, so this reader is
    the Perl's, followed line for line, and it has not been run against
    real data.
    """
    deck = Deck(path)
    parities, count = {}, 0
    for line in files.read_lines(path, strip=False):
        if not _LOPT.match(line):
            continue
        fields = line.split('\t')
        designation, j, energy, _, uncertainty = fields[:5]
        configuration, term = fields[8], fields[9].strip()
        if most is not None and _number(uncertainty) > most:
            continue
        parity = 'o' if designation.strip().endswith('*') else 'e'
        if not count:
            parities = {'e': 1, 'o': 2} if parity == 'e' \
                else {'e': 2, 'o': 1}
        count += 1
        ee = _number(energy) * KK
        vector = {'%s,%s' % (configuration.strip(), ' '.join(term.split())):
                  100}
        for one in _js(j, term):
            deck.add(Level(j=one, ee=ee, ec=ee, observed=True,
                           vector=dict(vector), kset=parities[parity]))
    return deck


def _number(text):
    try:
        return float(text)
    except ValueError:
        return 0.0


def _js(j, term):
    """The J values one LOPT row stands for.

    A row that names a J is one level.  A row with the J column empty is
    a term, and stands for every J the term has -- |L-S| to L+S.
    """
    if j.strip():
        return [j.strip()]
    found = re.match(r'^(\d+)([SPDFGHIKLMNOPQRSTUV])', term.strip())
    if not found:
        return []
    spin = (int(found.group(1)) - 1) / 2.0
    orbital = _L.index(found.group(2))
    out, value = [], abs(orbital - spin)
    while value <= orbital + spin + 1e-9:
        out.append(_j_text(value))
        value += 1.0
    return out


# --------------------------------------------------------------------------
# matching
# --------------------------------------------------------------------------

def distance(one, other):
    """How far apart two eigenvectors are, and whether that is certain.

    Returns `(distance, certain)`.  A component in one vector and not
    the other contributes its own square, so a level that shares nothing
    is as far away as the two vectors are long.
    """
    near = far = 0.0
    for key, out in other.items():
        if key not in one:
            near += out * out
            far += out * out
            continue
        inside = one[key]
        near += (inside - out) ** 2
        far += (inside + out) ** 2
        if abs(inside) > CERTAIN and abs(out) > CERTAIN:
            return 0.0, True
    for key, inside in one.items():
        if key not in other:
            near += inside * inside
            far += inside * inside
    return min(near, far), False


def identify(source, target):
    """Give every observed level of `source` a level of `target`.

    The candidates are tried in order of how near their fitted energy is
    to the source level's, so that where two levels of a J are equally
    good the closer one wins.  A target level that already carries an
    observed energy is passed over -- unless it is the very energy being
    assigned, which is how a deck already carrying some of these levels
    is brought up to date rather than refused.  Returns the list of
    source levels nothing could be found for.
    """
    lost = []
    for level in source.levels():
        if not level.observed:
            continue
        candidates = [one for one in target.group(level.kset, level.j)
                      if not one.taken
                      or abs(one.energy - level.ee) < 0.001]
        candidates.sort(key=lambda one: abs(one.ec - level.ec))
        best, least = None, None
        for one in candidates:
            apart, certain = distance(level.vector, one.vector)
            if certain:
                best = one
                break
            if least is None or apart < least:
                best, least = one, apart
        if best is None:
            lost.append(level)
        else:
            best.assigned = level.ee
            level.match = best
    return lost


# --------------------------------------------------------------------------
# what the run has to report
# --------------------------------------------------------------------------

def suspicious(source):
    """Matches worth a second look, as `print_suspicious` reports them.

    A match is quietly right where one of the two vectors is dominated
    by a single component, or where both vectors agree on their two
    largest.  Everything else is listed for the author to check.
    """
    out = []
    for level in source.levels():
        best = getattr(level, 'match', None)
        if not level.observed or best is None:
            continue
        mine, theirs = level.sorted_vector(), best.sorted_vector()
        if not mine or not theirs:
            continue
        if abs(mine[0][1]) > CERTAIN and abs(theirs[0][1]) > CERTAIN:
            continue
        if _same(mine, theirs, 2):
            continue
        out.append((level, best))
    return out


def _same(one, other, most):
    """Whether two component lists agree on their first `most` keys."""
    for i in range(most):
        left = one[i][0] if i < len(one) else None
        right = other[i][0] if i < len(other) else None
        if left != right:
            return False
    return True


def duplicates(target):
    """Target levels of one J given the same observed energy.

    Two levels cannot share an observed energy, so this is where a
    match has gone wrong; `print_duplicate` reports it the same way.
    """
    out = []
    for kset in sorted(target.sets):
        for j in sorted(target.sets[kset], key=_j_value):
            got = [one for one in target.group(kset, j)
                   if one.assigned is not None]
            got.sort(key=lambda one: (one.assigned, one.ec))
            for before, after in zip(got, got[1:]):
                if before.assigned == after.assigned and before.ec != after.ec:
                    out.append((before, after))
    return out


def known_configurations(target):
    """{parity set: [(configuration, lowest observed energy)]}.

    The configurations experiment has actually reached, lowest first --
    `print_confs`, which is what tells the author how far up the
    calculation the observations go.
    """
    out = {}
    for kset in sorted(target.sets):
        lowest = {}
        for j in target.sets[kset]:
            for one in target.group(kset, j):
                if one.assigned is None:
                    continue
                configuration = one.leading[1].partition(',')[0]
                if configuration not in lowest \
                        or one.assigned < lowest[configuration]:
                    lowest[configuration] = one.assigned
        out[kset] = sorted(lowest.items(), key=lambda kv: kv[1])
    return out


def vector_text(level, most=2):
    """A vector as the reports print it, in intensity per cent."""
    out = []
    for key, percent in level.sorted_vector()[:most]:
        configuration, _, term = key.partition(',')
        out.append('%4d%% %-7s%-8s' % (int(percent * percent / 100),
                                       configuration, term))
    return ''.join(out)
