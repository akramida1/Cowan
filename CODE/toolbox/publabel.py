"""Designations unique within each (parity, J).

The last of the three tools that turn a finished fit into a published
level list:

    levels.txt --12--> coupling_choice.txt --13--> levels_prepub.txt
                                                   --14--> levels_pub.txt

A published or database level needs a name no other level of the same
parity and J carries, because that triple is what a citation, a line
list and an ASD load all key on.  The leading component supplies the
name where a level is reasonably pure; where it is not, the name is
taken from further down the composition, and where that fails too, from
letter prefixes or from the energy.  That progression is not a loss:
the label of a strongly mixed level has little physical meaning in any
scheme, and the composition printed beside it is what carries the
physics.

WHAT IS PUBLISHED
-----------------

Truncation happens first.  Uniqueness is required only among the levels
that will actually be printed, and demanding it of a calculated tail
that no one will ever cite would spend letters and fallbacks on levels
nobody asked for.  The tail is cut per configuration set -- in practice
one set per parity, RCE being hard-limited to two (`CONFIG(IC,2)`,
`FOR/RCE20KD.F:170`) -- and the default cut is the position of the
highest observed level in the set, which is where the calculation stops
being anchored to measurement.

THE SEARCH
----------

Within one (parity, J) group each level offers its components at or
above `--min-name-percent` as candidates, and an assignment is wanted
that gives every level a different name.  This is a bipartite matching
and it is solved as one: the group is split into clusters that share no
candidate name -- nearly every level is a cluster of its own, since
nothing else in its group is named after the same basis state -- and
each cluster is searched exhaustively with a node budget.

The objective is the total chosen percentage with a strong preference
for the leading component (`LEADING`), under a much larger preference
for naming a level at all (`NAMED`).  There is no `>50%` pre-pass: a
level with a large leading component keeps it because the search
prefers it, which needs no threshold and cannot miss the case where two
levels of one group lead on the same basis state and the second-best
choice belongs to whichever of them loses less by it.

WHERE THE SEARCH CANNOT REACH
-----------------------------

A level the search leaves unnamed -- because every candidate it has is
taken, or because it has no candidate at all -- falls back to its
leading component with a letter prefix, `a `, `b `, ... on the term.

The letters are counted per (parity, configuration, term) across all J,
in increasing energy, and not within one (parity, J) group.  Within a
group the letter a level got would depend on how many same-named levels
happened to lie below it *in that J*, which changes from J to J as
levels cross, so one multiplet would come out `a 4F` at J=3/2 and
`b 4F` at J=5/2.  Counted over the family the letter is a property of
the level, and uniqueness inside a (parity, J) group holds either way,
since two levels of one group that want one name are two members of one
family and take two letters.

If more than `--letters-max` of the published levels need letters, or
more than `--leading-min` of them have a leading component under 50%,
then no naming scheme is going to mean anything and the tool stops
pretending: every level is named by the integer part of its energy,
with a blank term and a blank percentage, a sequence number having no
percentage to quote.  Should two integer energies of one parity collide,
a sequential index per parity is used instead.
"""

import string

from dataclasses import dataclass

from .purity import PurityError


#: The columns of levels_pub.txt before the composition groups.  `conf`
#: / `term` / `perc` are the chosen name and the weight it carries.
NAME = ['par', 'conf', 'term', 'perc']

#: What follows the groups: the energies, J and g `prepub` carries.
TAIL = ['Ee', 'Ec', 'dE_O-C', 'J', 'Lande_g']

#: cm-1 per unit of the energies as RCE writes them and as every file
#: up to `levels_prepub.txt` carries them.  This is the last file of the
#: three and the one a database loads, so it is the one that has to be
#: in the unit a spectroscopist quotes: the energies are multiplied up
#: on the way out, here and nowhere earlier (`iden.py:436` does the same
#: for the same reason).  1000 cm-1 is the unit in every deck of this
#: package -- `spectrum.py`, `radiate.py` and `explevs.KK` all say so --
#: so this is a conversion and not a guess.
CM1 = 1000.0

PUB = 'levels_pub.txt'

#: Worth of naming a level at all, against the worth of naming it well.
#: Larger than any attainable total of the rest, so an assignment that
#: names one more level always wins.
NAMED = 100000

#: Worth of naming a level after its own leading component.  Every level
#: is offered the same bonus, so where two of them want one name the
#: bonus cancels and the percentages decide which keeps it -- which is
#: the intended reading of "a strong preference for the leading
#: component" and not a veto on ever taking a lower one.
LEADING = 200

#: Nodes one cluster of the search may spend before it settles for the
#: best assignment it has found.  Reached only by a cluster of many
#: levels all competing for the same few names, where the difference
#: between the exhaustive answer and a very good one is a percent or
#: two of one label.
BUDGET = 200000


class PubError(PurityError):
    """levels_prepub.txt will not support the request."""


# --------------------------------------------------------------------------
# what gets published
# --------------------------------------------------------------------------

@dataclass
class Name:
    """A published level: the level, the name chosen, and where from."""

    level: object                   # prepub.Level
    conf: str = ''
    term: str = ''
    percent: object = None          # int, or None for an energy label
    index: int = -1                 # 0-based place in the composition
    letter: str = ''                # '', 'a', 'b', ...
    energy: bool = False            # named by its energy, not a component

    @property
    def label(self):
        """The name as it is printed, letter and all."""
        term = ('%s %s' % (self.letter, self.term)).strip() \
            if self.letter else self.term
        return (self.conf, term)

    @property
    def leading(self):
        return self.index == 0


def sets(levels):
    """The configuration sets, as the parities in the file's own order.

    A set is a parity here.  RCE takes two configuration sets and no
    more, so in practice they are one even and one odd, and
    levels_prepub.txt carries the parity but not the set number.  The
    order is the order the file meets them, which is the order of IN36's
    sets, which is the author's.
    """
    out = []
    for level in levels:
        if level.parity not in out:
            out.append(level.parity)
    return out


def cut(levels, nlev=None):
    """The levels that will be published, and how many from each set.

    `nlev` is one count per set, in the order `sets` gives them.  The
    default for a set is the position of its highest observed level.  A
    set with no observed level has no such position, and rather than
    publish all of it or none of it the run stops and names the set --
    only the author knows how far an unanchored calculation is worth
    printing.

    Returns `(published, counts)`, `counts` being {parity: n}.
    """
    order = sets(levels)
    if nlev is not None and len(nlev) != len(order):
        raise PubError(
            '--nlev takes one count per configuration set and this run has '
            '%d (%s), not %d'
            % (len(order), ', '.join(order), len(nlev)))

    counts = {}
    for at, parity in enumerate(order):
        if nlev is not None:
            counts[parity] = nlev[at]
            continue
        highest = 0
        for position, level in enumerate(
                [one for one in levels if one.parity == parity], 1):
            if not level.unfitted:
                highest = position
        if not highest:
            raise PubError(
                'the %s set has no observed level, so there is no highest '
                'observed level to cut at.  Say how many of its levels are '
                'to be published with --nlev (one count per set, in the '
                'order %s).' % (parity, '/'.join(order)))
        counts[parity] = highest

    out, so_far = [], {}
    for level in levels:
        so_far[level.parity] = so_far.get(level.parity, 0) + 1
        if so_far[level.parity] <= counts[level.parity]:
            out.append(level)
    return out, counts


# --------------------------------------------------------------------------
# the search
# --------------------------------------------------------------------------

def candidates(level, min_percent):
    """[(score, name, index)] for one level, best first."""
    out = []
    for index, component in enumerate(level.components):
        if component.percent < min_percent:
            continue
        score = component.percent + (LEADING if index == 0 else 0)
        out.append((score, (component.configuration, component.term), index))
    out.sort(key=lambda one: (-one[0], one[2]))
    return out


def _clusters(options):
    """Split the levels of a group by which of them share a name.

    Two levels that cannot want the same name cannot constrain each
    other, so they are searched apart.  In a real spectrum almost every
    level is a cluster of one and the search never runs at all.
    """
    parent = list(range(len(options)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a

    owner = {}
    for at, offered in enumerate(options):
        for _, name, _ in offered:
            if name in owner:
                union(owner[name], at)
            else:
                owner[name] = at
    out = {}
    for at in range(len(options)):
        out.setdefault(find(at), []).append(at)
    return list(out.values())


def _search(offered, budget=BUDGET):
    """Best distinct-name assignment for one cluster.

    Returns a list, one entry per level, each `None` or the level's
    `(score, name, index)`.  Most-constrained-first with an optimistic
    bound; a greedy pass seeds it so that the budget running out still
    leaves a usable answer rather than none.
    """
    order = sorted(range(len(offered)), key=lambda at: len(offered[at]))
    # The most any tail of `order` could still be worth, for the bound.
    best_left = [0] * (len(order) + 1)
    for at in range(len(order) - 1, -1, -1):
        one = offered[order[at]]
        best_left[at] = best_left[at + 1] + (NAMED + one[0][0] if one else 0)

    best = [None] * len(offered)
    taken = set()
    for at in order:                                    # the greedy seed
        for choice in offered[at]:
            if choice[1] not in taken:
                best[at] = choice
                taken.add(choice[1])
                break
    top = [sum(NAMED + one[0] for one in best if one)]
    spent = [0]
    chosen = [None] * len(offered)

    def step(depth, total, used):
        if spent[0] > budget:
            return
        spent[0] += 1
        if total + best_left[depth] <= top[0]:
            return
        if depth == len(order):
            top[0] = total
            best[:] = chosen
            return
        at = order[depth]
        for choice in offered[at]:
            if choice[1] in used:
                continue
            used.add(choice[1])
            chosen[at] = choice
            step(depth + 1, total + NAMED + choice[0], used)
            used.discard(choice[1])
        # Leaving a level unnamed is a last resort, and is tried last so
        # that the bound has a real total to prune against first.
        chosen[at] = None
        step(depth + 1, total, used)

    step(0, 0, set())
    return best


def choose(levels, min_percent=20.0, budget=BUDGET):
    """Name every level, one (parity, J) group at a time.

    Returns [(level, choice)] in the order given, `choice` being the
    level's `(score, name, index)` or `None` where the group had no name
    left for it.
    """
    groups = {}
    for at, level in enumerate(levels):
        groups.setdefault((level.parity, level.j), []).append(at)

    out = [None] * len(levels)
    for members in groups.values():
        offered = [candidates(levels[at], min_percent) for at in members]
        for cluster in _clusters(offered):
            if len(cluster) == 1:
                one = offered[cluster[0]]
                out[members[cluster[0]]] = one[0] if one else None
                continue
            answer = _search([offered[at] for at in cluster], budget=budget)
            for at, choice in zip(cluster, answer):
                out[members[at]] = choice
    return [(level, out[at]) for at, level in enumerate(levels)]


# --------------------------------------------------------------------------
# the letters, and the energy fallback
# --------------------------------------------------------------------------

def _letter(n):
    """0 -> 'a', 25 -> 'z', 26 -> 'aa'.  A spectrum needing 27 of these
    has already tripped `--letters-max`; the sequence goes on anyway
    rather than run out silently."""
    out = ''
    while True:
        out = string.ascii_lowercase[n % 26] + out
        n = n // 26 - 1
        if n < 0:
            return out


def name_all(chosen):
    """Turn the search's answer into `Name`s, letters and all.

    `chosen` is what `choose` returned.  Returns `(names, lettered,
    nameless)`: the names in the order given, the ones that needed a
    letter, and the levels that had no component to be named after at
    all and fell straight to their energy.
    """
    names, contested, nameless = [], set(), []
    for level, choice in chosen:
        if choice is not None:
            _, (conf, term), index = choice
            names.append(Name(level=level, conf=conf, term=term, index=index,
                              percent=level.components[index].percent))
            continue
        if not level.components:
            # Nothing survived phase 13's cut, so there is no basis
            # state to name it after and no letter can help.
            one = Name(level=level)
            names.append(one)
            nameless.append(one)
            continue
        # The search had no name left for it, so it takes its leading
        # component anyway and the whole contest is settled by letters.
        first = level.components[0]
        names.append(Name(level=level, conf=first.configuration,
                          term=first.term, index=0, percent=first.percent))
        contested.add((level.parity, level.j,
                       first.configuration, first.term))

    # Both sides of a contest take a letter, not just the one that lost
    # it: `4G*` beside `a 4G*` in one (parity, J) would read as two
    # different things, when they are two levels of one name that had to
    # be told apart.
    lettered = [one for one in names
                if (one.level.parity, one.level.j, one.conf, one.term)
                in contested]

    # The letters count over the (parity, configuration, term) family
    # across all J, in increasing energy -- which is the order of the
    # list, levels_prepub being written in Ec order within each parity.
    # Not within one (parity, J): there the letter a level got would
    # depend on how many same-named levels happened to lie below it in
    # that J, which changes from J to J as levels cross.
    used = {}
    for one in lettered:
        key = (one.level.parity, one.conf, one.term)
        one.letter = _letter(used.get(key, 0))
        used[key] = used.get(key, 0) + 1
    for one in nameless:
        _as_energy(one)
    return names, lettered, nameless


def _as_energy(one):
    energy = one.level.e_calc if one.level.unfitted else one.level.e_obs
    one.conf = '%d' % int(energy * CM1)
    one.term = ''
    one.percent = None
    one.index = -1
    one.letter = ''
    one.energy = True


def by_energy(names):
    """Name every level by its energy: the global fallback.

    The integer part of the observed energy where there is one and of
    the computed energy otherwise, with a blank term and a blank
    percentage.  Should two of one parity land on the same integer, the
    energies cannot separate the levels either and a sequential index
    per parity is used for the whole file -- one file, one convention,
    rather than half of it numbered.

    The integer part is of the energy in cm-1 (`CM1`), the unit this
    file is written in, so it is the whole number a spectroscopist would
    quote and it separates levels at the 1 cm-1 the list resolves.  Were
    it taken in the 10^3 cm-1 the earlier files carry, a thousand levels
    spread over 200 kK would collide almost everywhere and the
    sequential index would be the usual outcome rather than the rare
    one.
    """
    for one in names:
        _as_energy(one)
    seen, collides = set(), False
    for one in names:
        key = (one.level.parity, one.conf)
        if key in seen:
            collides = True
            break
        seen.add(key)
    if collides:
        at = {}
        for one in names:
            at[one.level.parity] = at.get(one.level.parity, 0) + 1
            one.conf = '%d' % at[one.level.parity]
    return collides


def fallback(names, lettered, letters_max=50.0, leading_min=50.0):
    """Whether names have stopped meaning anything.  (reason, share)

    Two ways for that to happen, and either is enough: most of the
    spectrum needing a letter, or most of it having no component large
    enough to be called dominant.  Both are shares of the published
    levels, and both default to half of them.
    """
    total = len(names) or 1
    share = 100.0 * len(lettered) / total
    if share > letters_max:
        return 'letters', share
    weak = sum(1 for one in names
               if not one.level.components
               or one.level.components[0].percent < 50)
    share = 100.0 * weak / total
    if share > leading_min:
        return 'mixing', share
    return None, 0.0


def outside(names, n_comp):
    """The levels named after a component the file does not print."""
    return [one for one in names if one.index >= n_comp]


# --------------------------------------------------------------------------
# levels_pub.txt
# --------------------------------------------------------------------------

def header(n_comp=3):
    return (NAME
            + [what % n for n in range(1, n_comp + 1)
               for what in ('perc%d', 'conf%d', 'term%d')]
            + TAIL)


def row_of(one, n_comp=3):
    """One line of levels_pub.txt, as column values.

    The group 1 columns are blank on every row named after its own
    leading component: that is the ASD rule -- a level says once where
    its name came from -- and it is what makes the file readable, the
    filled group 1 being exactly the rows whose name came from further
    down.  Reading rule: the leading component is group 1 where it is
    filled and `conf`/`term`/`perc` where it is not.
    """
    level = one.level
    conf, term = one.label
    out = [level.parity, conf, term,
           '' if one.percent is None else '%d' % one.percent]
    for at in range(n_comp):
        if at == 0 and one.leading:
            out.extend(['', '', ''])
            continue
        if at >= len(level.components):
            out.extend(['', '', ''])
            continue
        component = level.components[at]
        out.extend(['%d' % component.percent,
                    component.configuration, component.term])
    # In cm-1, not in the 10^3 cm-1 the rest of the pipeline carries:
    # this is the file that gets published.  Three decimals is what the
    # kK figures hold (6 decimals on Ee, 4 on Ec), so nothing is lost
    # and nothing invented.
    out.extend(['' if level.unfitted else '%.3f' % (level.e_obs * CM1),
                '%.3f' % (level.e_calc * CM1),
                '' if level.unfitted
                else '%.3f' % ((level.e_obs - level.e_calc) * CM1),
                '%4.1f' % level.j,
                '' if level.g is None else '%6.3f' % level.g])
    return out


def write_pub(path, names, n_comp=3):
    with open(path, 'w') as out:
        out.write('\t'.join(header(n_comp)) + '\n')
        for one in names:
            out.write('\t'.join(row_of(one, n_comp)) + '\n')


def duplicates(names):
    """Names that are not unique within their (parity, J).  Always [].

    The whole file is worth nothing if this is not empty, so it is
    checked rather than trusted, on the way out, every run.
    """
    seen, out = {}, []
    for one in names:
        key = (one.level.parity, one.level.j) + one.label
        if key in seen:
            out.append((seen[key], one))
        else:
            seen[key] = one
    return out
