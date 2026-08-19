"""Putting a level list together, and merging two runs into one.

Three jobs, in the order they have to happen:

  * **Index the basis.**  RCE names an eigenvector component by its
    configuration label and its state label -- `2p2  (3P) 4P` -- and RCG
    names the same thing by a number.  `BasisIndex` ties the two
    together, which is what makes everything below possible.
  * **Identify.**  RCE re-sorts the levels it fitted, so the n-th level
    of a J block in RCEOUT is generally *not* the n-th eigenvalue in
    OUTG11.  `identify()` pairs them by LS eigenvector composition, the
    method in `conv_cowan.pl:1182-1312`, and falls back to position only
    above the highest observed level, where RCE did not re-sort.
  * **Merge.**  A second run in a sub-folder supplies a coupling the
    first does not.  There are two kinds of second run and they must be
    matched differently -- see `match_mode()`.

Nothing here prints anything.  What a level is *called* is `asd.py`'s
business and what goes in a report is `printout.py`'s.
"""

from dataclasses import dataclass, field

from . import outg11, rce


class LevelsError(Exception):
    """Two lists of levels that cannot be put together."""


# --------------------------------------------------------------------------
# the basis of one J block
# --------------------------------------------------------------------------

@dataclass
class BasisIndex:
    """The LS basis states of one J block, by label and by number.

    A label is `(configuration, state)` exactly as both programs write
    it: RCG puts the state after an unmatched '(' in its row
    (`Row.bracketed`) and RCE puts the same text back with the '(' in
    front of it, so one lstrip is the whole translation.
    """

    numbers: dict = field(default_factory=dict)     # (cfg, state) -> number
    rows: list = field(default_factory=list)        # the section's rows
    _top: dict = field(default_factory=dict, repr=False)

    @classmethod
    def of(cls, section):
        numbers = {}
        for row in section.rows:
            key = (row.label, row.bracketed)
            # A repeated label would make the identification ambiguous;
            # RCG does not produce one, and if it ever did, guessing
            # would be worse than stopping.
            if key in numbers:
                raise LevelsError(
                    'two basis states of one J block carry the label '
                    '%s %s' % key)
            numbers[key] = row.number
        return cls(numbers=numbers, rows=list(section.rows))

    def number_of(self, component):
        """The basis-state number of an RCE eigenvector component."""
        return self.numbers.get((component.configuration,
                                 component.label.lstrip('(').rstrip()))

    def amplitude(self, state, level):
        """The amplitude of basis state `state` in eigenvector `level`.

        Both are 1-based, as RCG numbers them.  RCG prints amplitudes
        rather than percentages, so this is scaled by 100 to compare
        with what RCE prints.
        """
        row = self.rows[state - 1]
        if level > len(row.amplitudes):
            return 0.0
        return 100.0 * row.amplitudes[level - 1]

    def components_of(self, level, most=10):
        """An RCG eigenvector's largest components: [(amplitude, state)].

        Cached: identification asks for this once per candidate and a
        Pr III block has three hundred of them.
        """
        cached = self._top.get(level)
        if cached is None:
            cached = [(self.amplitude(row.number, level), row.number)
                      for row in self.rows]
            cached.sort(key=lambda pair: -abs(pair[0]))
            cached = cached[:10]
            self._top[level] = cached
        return cached[:most]


def by_block(run, couplings=('LS',)):
    """{(parity set, J): JBlock} for a whole OUTG11."""
    out = {}
    for block in run.blocks(list(couplings)):
        out[(block.kset, block.j)] = block
    return out


# --------------------------------------------------------------------------
# identification
# --------------------------------------------------------------------------

@dataclass
class Identified:
    """One RCE level, tied to the RCG eigenvector it came from."""

    level: object                   # rce.Level
    kset: int
    j: float
    rcg: int                        # eigenvector number in its J block
    distance: float = None          # None where position was used
    positional: bool = False


def eigenvalue(blocks, one):
    """The RCG eigenvalue an identified level came from.

    This is the ordering key of every printed level list, on the
    author's instruction of 2026-08-15: RCE prints its own list in the
    order of its Ecalc, and RCG's eigenvalue is the same quantity
    recomputed from the rounded parameters RCE punched back into ING11.
    OUTG11 prints it to three decimals, so the fitted energy is returned
    beside it to break a tie deterministically -- `printout.pl` left
    such ties in the order Perl's hash iteration happened to produce.
    """
    block = blocks.get((one.kset, one.j))
    if block is not None and one.rcg <= len(block.eigenvalues):
        return (block.eigenvalues[one.rcg - 1], one.level.e_calc)
    return (one.level.e_calc, one.level.e_calc)


def unfitted(blocks, kset):
    """Stand-in levels for a run with no RCE output: RCG's own eigenvalues.

    Everything downstream works on `Identified`, so a run that was never
    fitted is given one per eigenvector, with the eigenvalue standing in
    for both the observed and the fitted energy and the Lande g factor
    taken from the block where RCG printed one.
    """
    out = []
    for (k, j), block in sorted(blocks.items()):
        if k != kset:
            continue
        for number, energy in enumerate(block.eigenvalues, 1):
            g = (block.g_values[number - 1]
                 if number <= len(block.g_values) else None)
            level = rce.Level(e_obs=energy, e_calc=energy, j=j, flag='*',
                              configuration='', g=g, components=(),
                              kset=kset, index=number)
            out.append(Identified(level=level, kset=kset, j=j, rcg=number,
                                  positional=True))
    return out


def identify(cset, blocks, tolerance=0.005, emax=None, skipped=None):
    """Pair every level of one RCE parity set with an RCG eigenvector.

    RCE sorts the levels it has observations for; above the highest of
    those it leaves RCG's own order alone.  So the levels below that
    energy are matched on composition and the rest by position, which
    is exactly what `conv_cowan.pl:1212-1305` does -- and the boundary
    is the highest *fitted* energy among the observed levels, not the
    highest observed one.

    `emax` is the deck's own E0MAX (`ing11.Ing11.emax`).  A level above
    it that cannot be matched is left out of the result and appended to
    `skipped` instead of stopping the run: RCG prints nothing above
    E0MAX, so RCE may well list a level whose eigenvector is not in the
    file, and that is the deck saying so rather than a failure.  This is
    the condition at `conv_cowan.pl:1299`, where an unmatched level dies
    unless E0MAX is set and the level is above it.  With `emax` left
    None -- what every other caller passes -- nothing changes.
    """
    out = []
    for j, levels in sorted(cset.by_j().items()):
        block = blocks.get((cset.index, j))
        if block is None:
            raise LevelsError(
                'OUTG11 has no J=%s block in configuration set %d, but '
                'RCEOUT lists %d levels for it'
                % (_j(j), cset.index, len(levels)))
        index = BasisIndex.of(block.section('LS'))
        ceiling = _ceiling(levels, tolerance)
        taken = set()
        for position, level in enumerate(levels, 1):
            if ceiling is None or level.e_calc >= ceiling:
                if position in taken:
                    raise LevelsError(
                        'RCE level %d (set %d, J=%s, E=%.3f) falls on RCG '
                        'eigenvector %d, which is already taken'
                        % (position, cset.index, _j(j), level.e_calc,
                           position))
                taken.add(position)
                out.append(Identified(level=level, kset=cset.index, j=j,
                                      rcg=position, positional=True))
                continue
            best, distance = _best_match(level, block, index, taken)
            if best is None and emax and level.e_calc > emax:
                if skipped is not None:
                    skipped.append(level)
                continue
            if best is None:
                raise LevelsError(
                    'no RCG eigenvector matches the RCE level at E=%.3f '
                    '(set %d, J=%s).  Its components are %s'
                    % (level.e_calc, cset.index, _j(j),
                       ', '.join('%s %s' % (c.configuration, c.label)
                                 for c in level.components)))
            taken.add(best)
            out.append(Identified(level=level, kset=cset.index, j=j,
                                  rcg=best, distance=distance))
    return out


def _ceiling(levels, tolerance):
    """The highest fitted energy among the levels RCE had an observation
    for; below it RCE's order is its own, above it RCG's."""
    observed = [l.e_calc for l in levels if l.observed]
    return max(observed) + tolerance if observed else None


def _best_match(level, block, index, taken):
    """The RCG eigenvector closest to one RCE level's composition.

    The candidates are the eigenvectors that carry at least one of the
    level's own components; they are tried nearest-in-energy first, and
    the distance is the sum of squared differences between the ten
    printed percentages and RCG's own, computed with both relative signs
    because an eigenvector is defined only up to one.  Components RCG
    has and RCE did not print add their whole square, up to four of
    them, which is what keeps a small overlap from winning.
    """
    wanted = []
    for component in level.components:
        state = index.number_of(component)
        if state is not None:
            wanted.append((state, float(component.percent)))
    if not wanted:
        return None, None

    candidates = set()
    for state, _ in wanted:
        row = index.rows[state - 1]
        for n in range(1, len(row.amplitudes) + 1):
            candidates.add(n)
    order = sorted(candidates,
                   key=lambda n: abs(block.eigenvalues[n - 1] - level.e_calc)
                   if n <= len(block.eigenvalues) else float('inf'))

    best, least = None, None
    for n in order:
        if n in taken:
            continue
        same = other = 0.0
        for state, percent in wanted:
            a = index.amplitude(state, n)
            same += (a - percent) ** 2
            other += (a + percent) ** 2
        seen = {state for state, _ in wanted}
        extra = 0
        for amplitude, state in index.components_of(n):
            if state in seen:
                continue
            same += amplitude * amplitude
            other += amplitude * amplitude
            extra += 1
            if extra >= 4:
                break
        distance = min(same, other)
        if least is None or distance < least:
            best, least = n, distance
    return best, least


def _j(j):
    return '%d/2' % round(2 * j) if round(2 * j) % 2 else '%g' % j


# --------------------------------------------------------------------------
# merging a second run
# --------------------------------------------------------------------------

EIGENVECTOR = 'eigenvector'
ENERGY_ORDER = 'energy-order'


def match_mode(primary, secondary):
    """How two decks' levels may be matched: by composition or by order.

    Same configuration cards in the same slot order means the basis
    states are the same objects, so the two runs' LS eigenvectors may be
    compared directly.  The same cards with the slots *permuted* is a
    reordered deck: the basis differs, nothing can be compared, and the
    only thing left is the energy order within a (parity, J) group.
    Different cards altogether is not a pair of runs of one case.
    """
    a = [(c.kset, c.label, tuple(c.slots)) for c in primary]
    b = [(c.kset, c.label, tuple(c.slots)) for c in secondary]
    if a == b:
        return EIGENVECTOR
    loose_a = sorted((k, l, tuple(sorted(s))) for k, l, s in a)
    loose_b = sorted((k, l, tuple(sorted(s))) for k, l, s in b)
    if loose_a == loose_b:
        return ENERGY_ORDER
    only_a = sorted(set(l for _, l, _ in a) - set(l for _, l, _ in b))
    only_b = sorted(set(l for _, l, _ in b) - set(l for _, l, _ in a))
    raise LevelsError(
        'the two decks are not the same case.\n'
        '  only in the first:  %s\n  only in the second: %s'
        % (', '.join(only_a) or '(none)', ', '.join(only_b) or '(none)'))


@dataclass
class Merged:
    """One level of the primary run, with its twin in the secondary."""

    primary: object                 # Identified
    secondary: object               # Identified, or None where none matched
    delta: float = None             # E(secondary) - E(primary), energy order

    @property
    def matched(self):
        return self.secondary is not None


def merge_by_energy_order(primary, secondary):
    """Pair the n-th level of each (parity, J) with the n-th of the other.

    This is the only thing available when the two decks sum their
    subshells in different orders: the basis states differ, so no
    eigenvector can be compared.  If the reordering were exact the
    energies would be identical; in practice they differ a little, and
    that difference is reported rather than hidden -- see
    `energy_order_report()`.
    """
    groups = _grouped(secondary)
    out = []
    for key, ours in sorted(_grouped(primary).items()):
        theirs = groups.get(key, [])
        for n, one in enumerate(ours):
            other = theirs[n] if n < len(theirs) else None
            delta = (other.level.e_calc - one.level.e_calc) if other else None
            out.append(Merged(primary=one, secondary=other, delta=delta))
    return out


def merge_by_eigenvector(primary, secondary, blocks, other_blocks):
    """Pair two runs that share a basis, on LS eigenvector composition.

    Both runs carry an LS section for every block, so a level of one is
    matched to the level of the other whose LS eigenvector it is -- via
    the RCG eigenvector each was already identified with, which is the
    same object in both runs when the decks agree.

    Agreeing configuration cards are not enough for that: a term cutoff
    (NLT11, or columns 83-85 of a configuration card) shrinks the basis
    without changing a single card, so the block sizes are checked here
    and a difference stops the merge rather than silently pairing two
    different eigenvectors.
    """
    _check_basis(blocks, other_blocks)
    by_key = {}
    for one in secondary:
        by_key[(one.kset, one.j, one.rcg)] = one
    out = []
    for one in primary:
        other = by_key.get((one.kset, one.j, one.rcg))
        delta = (other.level.e_calc - one.level.e_calc) if other else None
        out.append(Merged(primary=one, secondary=other, delta=delta))
    return out


def _check_basis(blocks, other_blocks):
    """Both runs must have the same number of basis states per block."""
    bad = []
    for key, block in sorted(blocks.items()):
        other = other_blocks.get(key)
        if other is None:
            bad.append('set %d J=%s is missing from the second run'
                       % (key[0], _j(key[1])))
        elif block.nlevels != other.nlevels:
            bad.append('set %d J=%s has %d basis states in the first run and '
                       '%d in the second'
                       % (key[0], _j(key[1]), block.nlevels, other.nlevels))
    if bad:
        raise LevelsError(
            'the two runs do not share a basis, so their eigenvectors '
            'cannot be compared:\n  ' + '\n  '.join(bad[:8])
            + ('\n  ... and %d more' % (len(bad) - 8) if len(bad) > 8 else '')
            + '\nA term cutoff is the usual cause: check NLT11 (columns '
              '26-27 of the ING11 main card) and columns 83-85 of the '
              'configuration cards in both decks.')


def _grouped(identified):
    """{(parity set, J): [Identified]} in energy order."""
    out = {}
    for one in identified:
        out.setdefault((one.kset, one.j), []).append(one)
    for group in out.values():
        group.sort(key=lambda x: x.level.e_calc)
    return out


def energy_order_report(merged):
    """The largest |dE| per (parity set, J), worst first.

    A reordered run should reproduce the parent's energies; a group
    where it does not is either a genuinely broken reordering or a
    mismatched pair, and either way it has to be visible.
    """
    worst = {}
    for one in merged:
        if one.delta is None:
            continue
        key = (one.primary.kset, one.primary.j)
        if key not in worst or abs(one.delta) > abs(worst[key].delta):
            worst[key] = one
    rows = [(abs(m.delta), key, m) for key, m in worst.items()]
    rows.sort(reverse=True)
    return [(key, m) for _, key, m in rows]


def unmatched(merged):
    """Every level of the primary that found no twin."""
    return [m for m in merged if not m.matched]
