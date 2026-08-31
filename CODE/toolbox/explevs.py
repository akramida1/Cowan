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

from dataclasses import dataclass, field

from . import files, rce

# A component below this is not worth carrying: `subst_exp_levs.pl:110`
# keeps only those at or past 10, and stops at the first one that is not.
FLOOR = 10

# Both vectors holding this much of one basis state settle it
# (`subst_exp_levs.pl:308`).  70 in amplitude is 49% in intensity, and
# two components cannot both be that large.
CERTAIN = 70

# The same thing in intensity per cent, which is where a source that
# is not an RCE deck is compared: 70 in amplitude is 49 in intensity.
CERTAIN_INTENSITY = 49

# RCE prints five components on a level line and `subst_exp_levs.pl:74-86`
# reads all five.
COMPONENTS = 5

#: How a source level may be tied to a target level.
MATCH_BY = ('heuristic', 'ecalc', 'eexp')

# --------------------------------------------------------------------------
# the match key
#
# A deck-to-deck run compares the keys RCE itself writes, `conf,term`,
# and nothing below is used.  A source that is not a deck names its
# components with whatever it has -- a configuration, the parent term
# written inside it, the final term -- and any of the three may be
# missing, so the key it is matched on is the three fields separately
# with `ANY` where the source said nothing.  `project` then adds the
# target's components together over the fields the source left open,
# which is what lets one source name a component to the parent and
# another only to the configuration and both be compared fairly.
# --------------------------------------------------------------------------

#: A field the source (or the target's own label) does not pin down.
ANY = '*'

#: What separates the three fields of a match key.  Not a comma: the
#: deck's own keys use one, and the two spaces must not be confusable.
KEY = '|'


def match_key(configuration, parent, term):
    """The key a component is matched on, `ANY` for what is not said."""
    return KEY.join(field or ANY for field in (configuration, parent, term))


def key_fields(key):
    """The three fields of a match key, or `None` for a deck's own key."""
    if KEY not in key:
        return None
    return tuple(key.split(KEY, 2))


def key_matches(one, other):
    """Whether two match keys can stand for the same basis state."""
    for mine, theirs in zip(one.split(KEY), other.split(KEY)):
        if mine != theirs and mine != ANY and theirs != ANY:
            return False
    return True


def _sharpness(key):
    """How many of a key's fields are pinned down; the most specific first."""
    return -sum(1 for field in key.split(KEY) if field != ANY)


def project(mine, theirs):
    """The target's vector added up over the fields the source left open.

    A source that names `4f2.(3F).6p` and no final term is saying
    *this parent*, and the target's share of it is every component of
    that configuration and parent whatever their final term.  So the
    target is projected onto the keys the source actually used, most
    specific key first so that a component claimed by a sharp key is
    not claimed again by a blunt one.
    """
    out, taken = {}, set()
    for key in sorted(mine, key=_sharpness):
        total = 0.0
        for other, value in theirs.items():
            if other in taken or not key_matches(key, other):
                continue
            total += value
            taken.add(other)
        out[key] = total
    return out

class ExpLevelsError(Exception):
    """Two decks that cannot be matched level for level."""


@dataclass
class Level:
    """One level of an RCE deck, as the matcher needs it."""

    j: str                          # as written: '3', '3/2'
    ee: float                       # observed energy, or the computed one
    ec: float                       # fitted energy
    observed: bool
    status: str = ' '               # column 14 as the deck spells it: blank,
                                    # '*', or a group letter -- see `rce`
    vector: dict = field(default_factory=dict)   # {'conf,term': percent}
    kset: int = 1
    position: int = 0               # its place within its J, 1-based
    line: int = -1                  # the line it was read from
    assigned: float = None          # the energy this level has been given
    assigned_status: str = ' '      # the status character that came with it
    match: object = None            # the level of the other deck it is tied to
    mvector: dict = None            # the vector to match on, where it is not
                                    # `vector` -- see `altlevsrc`
    intensity: bool = False         # whether `vector` is already in intensity
    pvector: dict = None            # a target level's vector as the source
                                    # that took it saw it -- see `project`

    @property
    def match_vector(self):
        """The vector the matcher compares, `vector` unless told otherwise.

        A source that is not an RCE deck quotes intensity per cent and
        names its components the way a table does, so the target's
        amplitudes have to be brought into that space before the two
        can be compared; `altlevsrc.normalize_target` puts the result
        here and leaves `vector` as RCE wrote it, which is what the
        reports print.
        """
        return self.vector if self.mvector is None else self.mvector

    @property
    def group(self):
        """The heteroscedastic group this level belongs to, or `''`."""
        return rce.group_letter(self.status)

    @property
    def absolute(self):
        """Whether the level's observed energy is on the absolute scale."""
        return self.observed and not self.group

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

    def sorted_match_vector(self):
        """The components the matcher sees, largest first.

        For a target level tied to a source that named its components
        loosely, that is the projection the match was decided on and
        not the finer vector the deck holds, so that the report prints
        the two sides in the one key space they were compared in.
        """
        vector = self.match_vector if self.pvector is None else self.pvector
        return sorted(vector.items(), key=lambda kv: -abs(kv[1]))


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
                           status=level.flag or ' ',
                           vector=_vector(level.components), kset=cset.index,
                           line=places.get((cset.index, level.index), -1)))
    return deck


# --------------------------------------------------------------------------
# matching
# --------------------------------------------------------------------------

def distance(one, other, certain=CERTAIN):
    """How far apart two eigenvectors are, and whether that is certain.

    Returns `(distance, certain)`.  A component in one vector and not
    the other contributes its own square, so a level that shares nothing
    is as far away as the two vectors are long.

    `certain` is the share of one basis state that settles a match, and
    it depends on which space the comparison is being made in: 70 in
    amplitude per cent, which is RCE's, and `CERTAIN_INTENSITY` in the
    intensity per cent an experimental table quotes.
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
        if abs(inside) > certain and abs(out) > certain:
            return 0.0, True
    for key, inside in one.items():
        if key not in other:
            near += inside * inside
            far += inside * inside
    return min(near, far), False


def distance_named(one, other, certain=CERTAIN_INTENSITY):
    """The same, for a source that names only some of its components.

    A published table quotes two or three components and stops, so the
    target has ten and the source has two, and `distance` -- which
    charges every component one side has and the other does not -- then
    prefers whichever target level is most thoroughly mixed, because
    that is the one whose components are all small.  It is the wrong
    question: what a partial source says is *these* components and this
    much of something else, so that is what is compared -- the named
    components one by one, and everything the source did not name as a
    single remainder on both sides.

    Where the source does account for the whole level the two agree, so
    this is what an intensity comparison uses throughout.
    """
    apart = mine_total = theirs_total = 0.0
    for key, mine in one.items():
        theirs = other.get(key, 0.0)
        apart += (mine - theirs) ** 2
        mine_total += mine
        theirs_total += theirs
        if mine > certain and theirs > certain:
            return 0.0, True
    apart += (max(0.0, 100.0 - mine_total)
              - max(0.0, 100.0 - theirs_total)) ** 2
    return apart, False


def identify(source, target, by='heuristic', certain=CERTAIN):
    """Give every observed level of `source` a level of `target`.

    Three ways, `MATCH_BY`:

    * **heuristic**, the default -- a source level that carries a
      composition is tied to the candidate whose composition is
      nearest, and one that does not is tied on energy alone.  This is
      what a deck-to-deck run has always done, and it is what a source
      file with only some of its levels designated needs.
    * **ecalc** -- energy alone, against the target's fitted energies,
      ignoring any composition.
    * **eexp** -- energy alone, against the target's *observed*
      energies, and only levels that already carry one are candidates.
      That is how a deck built early in an analysis has its observed
      energies brought up to date from a later level optimization: the
      adjustments are small, so nearest wins, and where two source
      levels of one (parity, J) land on the same target level the
      duplicate report is what catches it.

    Candidates are tried in order of how near their energy is to the
    source level's, so that where two are equally good the closer one
    wins.  A target level that already carries an observed energy is
    passed over -- unless it is the very energy being assigned, which
    is how a deck already holding some of these levels is brought up to
    date rather than refused.  Returns the source levels nothing could
    be found for.
    """
    lost = []
    for level in source.levels():
        if not level.observed:
            continue
        group = target.group(level.kset, level.j)
        if by == 'eexp':
            candidates = [one for one in group
                          if one.observed and one.assigned is None]
            candidates.sort(key=lambda one: abs(one.ee - level.ee))
        else:
            candidates = [one for one in group
                          if not one.taken
                          or abs(one.energy - level.ee) < 0.001]
            candidates.sort(key=lambda one: abs(one.ec - level.ec))
        mine = level.match_vector
        seen = None
        if by != 'heuristic' or not mine:
            best = candidates[0] if candidates else None
        else:
            apply = distance_named if level.intensity else distance
            best, least = None, None
            for one in candidates:
                # A source that is not a deck names its components to
                # whatever depth it can, so the target is added up over
                # the fields it left open before the two are compared.
                theirs = (project(mine, one.match_vector) if level.intensity
                          else one.match_vector)
                apart, sure = apply(mine, theirs, certain)
                if sure:
                    best, seen = one, theirs
                    break
                if least is None or apart < least:
                    best, least, seen = one, apart, theirs
        if best is None:
            lost.append(level)
        else:
            best.assigned = level.ee
            # The status comes across with the energy.  A group level
            # whose letter were dropped here would be written into the
            # target as an ordinary absolute one, which is the silent
            # promotion `03_design.md` Sec. 7 asks `set_observed` to
            # make impossible.
            best.assigned_status = level.status if level.observed else ' '
            best.pvector = seen if level.intensity else None
            level.match = best
    return lost


# --------------------------------------------------------------------------
# what the run has to report
# --------------------------------------------------------------------------

def suspicious(source, certain=CERTAIN):
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
        mine = level.sorted_match_vector()
        theirs = best.sorted_match_vector()
        if not mine or not theirs:
            continue
        if abs(mine[0][1]) > certain and abs(theirs[0][1]) > certain:
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


def key_text(key):
    """(configuration, term) of a key, however it was written.

    A deck's own key is `conf,term`; a match key is the three fields of
    `match_key`, and its parent comes back inside the term the way RCE
    prints one, so that the two sides of a report read alike.
    """
    fields = key_fields(key)
    if fields is None:
        configuration, _, term = key.partition(',')
        return configuration, term
    configuration, parent, term = fields
    if parent != ANY:
        term = ('(%s) %s' % (parent, '' if term == ANY else term)).strip()
    elif term == ANY:
        term = ''
    return ('' if configuration == ANY else configuration), term


def vector_text(level, most=2):
    """A vector as the reports print it, in intensity per cent.

    RCE writes a signed amplitude per cent, so it is squared here; a
    source that is not an RCE deck quotes the intensity already and
    says so with `intensity`, and is printed as it came.
    """
    out = []
    for key, percent in level.sorted_vector()[:most]:
        configuration, term = key_text(key)
        out.append('%4d%% %-7s%-8s'
                   % (int(percent if level.intensity
                          else percent * percent / 100),
                      configuration, term))
    return ''.join(out)
