"""One coupling scheme per configuration, applied to every level.

Phase 12 measured how well each scheme describes each configuration and
the author marked a choice in `coupling_choice.txt`.  This is where the
choice is used, and the whole tool is one rule:

    for each level, take from each scheme's block only those components
    whose configuration was assigned that scheme; concatenate across
    schemes; sort descending by percentage; print.

So a level may read `44% 4f2.(3P*<0>).6p (0,3/2)*` followed by
`35% 4f.(2F*).5d2.(1G) 2D*` -- the first configuration written in JJ,
the second in LS, in one composition.  That is not sloppiness.  A change
of coupling recouples angular momenta *within* a configuration, so the
recoupling matrix is block-diagonal in the configuration and the weight
a configuration carries in a level is the same in every scheme
(`toolbox/purity.py`).  Only the spread of that weight over the
configuration's own basis states moves.  The percentages therefore still
add to the same total whichever scheme each configuration is quoted in,
and Sugar's Pr III tables are written exactly this way.

That invariance is checked rather than assumed: `disagreements()`
compares each configuration's summed weight between the schemes, level
by level, and reports where two schemes differ by more than the printing
itself can account for.

One thing the rule cannot repair.  The components come from
`printout.py`'s per-block cut, so a configuration's components may have
been cut out of the very block it was assigned while surviving in
another.  Run `printout --components 10 --min-percent 0.5` for this
pipeline and let `prepub` do the trimming; `shortfalls()` reports the
levels where the cut still bit, measured the only way that means
anything -- against what a single block already covered.
"""

from dataclasses import dataclass, field

from .purity import Component, PurityError


#: The columns of levels_prepub.txt before the components.  No `cpl`:
#: one line now carries every scheme, and the scheme a component is
#: written in is implied by its label.  No `e_c`: an unfitted level's
#: `Ee` is blank, which says the same thing without a second column.
#:
#: `suffix` is the compilation's own `+x`, published beside the energy
#: rather than appended to it: `Ee` stays a bare number, so `dE_O-C` and
#: every spreadsheet that reads the column keep working, and the suffix
#: is a field of its own that a reader can sort and filter on.  It is
#: written with a leading space -- `' +x'` -- because these files are
#: opened in Excel, where a cell beginning `+x` is parsed as a formula
#: referring to an undefined name.
FIXED = ['par', 'Ee', 'suffix', 'Ec', 'dE_O-C', 'J', 'Lande_g']


def suffix_text(letter):
    """`'x'` -> `' +x'`; `''` -> `''`.  The leading space is required."""
    return ' +%s' % letter if letter else ''

PREPUB = 'levels_prepub.txt'


# --------------------------------------------------------------------------
# the levels, gathered out of their blocks
# --------------------------------------------------------------------------

@dataclass
class Level:
    """One level, with its per-coupling blocks and the merged whole."""

    parity: str
    e_obs: float
    unfitted: bool
    e_calc: float
    j: float
    g: object
    conf0: str
    line: int                                       # of its first block
    group: str = ''                                 # the `+x` system, or ''
    blocks: dict = field(default_factory=dict)      # coupling -> Row
    components: list = field(default_factory=list)  # the merged composition

    @property
    def key(self):
        return (self.parity, self.e_calc, self.j)


def _key(row):
    return (row.parity, row.e_calc, row.j)


def collect(table):
    """Pair the rows of every coupling into one `Level` each.

    A levels.txt is the same levels over and over, one block per
    coupling, so they pair on (parity, Ec, J) -- the fields no change of
    coupling can move.  Two levels of one parity and J at the same
    computed energy would collide, so occurrences are counted and the
    n-th of one block pairs with the n-th of another; the blocks are
    written in eigenvalue order by one loop, so their n-th occurrences
    are the same level.

    Returns `(levels, missing)`, `missing` being {coupling: [Level]} for
    the levels a block does not carry -- which a merge from a second
    folder can produce, and which is reported, never dropped in silence.
    """
    levels, at, seen = [], {}, {}
    for row in table.rows:
        index = seen.get((row.coupling, _key(row)), 0)
        seen[(row.coupling, _key(row))] = index + 1
        slot = (_key(row), index)
        level = at.get(slot)
        if level is None:
            level = Level(parity=row.parity, e_obs=row.e_obs,
                          unfitted=row.unfitted, e_calc=row.e_calc,
                          j=row.j, g=row.g, conf0=row.conf0, line=row.line,
                          group=row.group)
            at[slot] = level
            levels.append(level)
        if row.coupling in level.blocks:
            raise PurityError(
                '%s line %d: the %s block gives the level (%s, Ec=%s, J=%s) '
                'twice' % (table.path, row.line, row.coupling,
                           row.parity, row.e_calc, row.j))
        level.blocks[row.coupling] = row
        # A level whose first block lost its optional fields still takes
        # them from a later one; every block prints the same energies,
        # the same g and the same conf0.
        if level.g is None:
            level.g = row.g
        if not level.conf0:
            level.conf0 = row.conf0

    missing = {}
    for coupling in table.couplings:
        absent = [one for one in levels if coupling not in one.blocks]
        if absent:
            missing[coupling] = absent
    return levels, missing


# --------------------------------------------------------------------------
# applying the choice
# --------------------------------------------------------------------------

def apply_choice(levels, choice, couplings, deck=None):
    """Build each level's composition out of the assigned blocks.

    `choice` is {configuration label: coupling} as `read_choice` gives
    it.  A component's configuration is its `conf` column stripped of
    the parenthesised intermediate terms.

    Returns `unassigned`, {label: count}, for the configurations that
    appear in a composition and that `coupling_choice.txt` does not
    name.  `config_purity.py` writes a line only for a configuration
    that dominates at least one level, so a configuration that is never
    more than a small admixture legitimately has no line; those take the
    first scheme of the choice file and are reported rather than
    refused, because dropping them would take their weight out of the
    level and the percentages would stop adding up.  A label that is in
    no scheme *and* in no deck is the other case, and that one stops the
    run by name: levels.txt and IN36 are then not from one run.
    """
    fallback = couplings[0]
    unassigned = {}
    known = None if deck is None else set(deck)
    for level in levels:
        merged = []
        for coupling in couplings:
            row = level.blocks.get(coupling)
            if row is None:
                continue
            for order_in_block, component in enumerate(row.components):
                label = component.bare
                wanted = choice.get(label)
                if wanted is None:
                    if known is not None and label not in known:
                        raise PurityError(
                            '%s line %d: the configuration %r of the '
                            'component %r is not one of the deck, so no '
                            'coupling scheme can be chosen for it; '
                            'levels.txt and IN36 are not from the same run.'
                            % (coupling, row.line, label,
                               component.configuration))
                    unassigned[label] = unassigned.get(label, 0) + 1
                    wanted = fallback
                if wanted == coupling:
                    merged.append((component, couplings.index(coupling),
                                   order_in_block))
        # Descending by percentage; ties keep the order of the schemes
        # and, within one scheme, the order its block printed them in,
        # so the file is reproducible rather than sorted by accident.
        merged.sort(key=lambda item: (-item[0].percent, item[1], item[2]))
        level.components = [component for component, _, _ in merged]
    return unassigned


def trim(levels, components=5, min_percent=1.5):
    """Cut every composition down to what is to be published."""
    for level in levels:
        kept = [c for c in level.components if c.percent >= min_percent]
        level.components = kept[:components]


# --------------------------------------------------------------------------
# the self-check
# --------------------------------------------------------------------------

@dataclass
class Disagreement:
    """One configuration whose weight in one level moved with the scheme."""

    level: Level
    label: str
    weights: dict                   # coupling -> summed percentage
    slack: float                    # what the printing can account for

    @property
    def spread(self):
        values = list(self.weights.values())
        return max(values) - min(values)


def disagreements(levels, couplings):
    """Where a configuration's weight in a level is not scheme-invariant.

    The weight is invariant in exact arithmetic; what is printed is not
    exactly the weight.  A block loses whatever fell below `printout`'s
    cut or past its component count, and every percentage is an integer
    rounded on its own, so a sum over k of them carries up to k/2 of
    rounding -- which is why a configuration spread over ten components
    can print as 102%.  Two blocks may therefore honestly differ by

        sum over the two blocks of (what the block did not print
                                    + half its count of components
                                      in this configuration)

    and no further.  A configuration that exceeds that is worth looking
    at, and is what this returns; nothing is corrected here.
    """
    out = []
    for level in levels:
        shortfall, weight, count = {}, {}, {}
        for coupling in couplings:
            row = level.blocks.get(coupling)
            if row is None:
                continue
            total = 0
            for component in row.components:
                per = weight.setdefault(component.bare, {})
                per[coupling] = per.get(coupling, 0) + component.percent
                many = count.setdefault(component.bare, {})
                many[coupling] = many.get(coupling, 0) + 1
                total += component.percent
            shortfall[coupling] = max(0, 100 - total)
        for label, per in weight.items():
            if len(per) < 2:
                # A configuration only one block printed says nothing:
                # the others cut it, they did not disagree about it.
                continue
            slack = sum(shortfall[c] + 0.5 * count[label][c] for c in per)
            values = list(per.values())
            if max(values) - min(values) > slack:
                out.append(Disagreement(level=level, label=label,
                                        weights=dict(per), slack=slack))
    return out


def shortfalls(levels, limit=1):
    """Levels the merge covers less well than a single block already did.

    Not a comparison with 100%.  A level of a 2500-level hole spectrum
    whose leading component is 11% has most of its weight below
    `printout`'s cut in *every* scheme, and the merge is not what lost
    it.  What the merge can lose is a configuration cut out of the very
    block it was assigned while it survived in another, and the symptom
    of that is a merged composition that covers less than one of the
    blocks it was built from -- the merge should normally cover more
    than any of them, taking each configuration from where it is best
    concentrated.

    Returns [(level, merged total, best block total)], measured before
    `trim` deliberately shortens anything.
    """
    out = []
    for level in levels:
        total = sum(c.percent for c in level.components)
        best = max([sum(c.percent for c in row.components)
                    for row in level.blocks.values()] or [0])
        if best - total > limit:
            out.append((level, total, best))
    return out


# --------------------------------------------------------------------------
# levels_prepub.txt
# --------------------------------------------------------------------------

def order(levels):
    """Parity blocks in the order levels.txt met them, then by Ec.

    Not sorted on the parity letter: the parities come out of RCE in the
    order of IN36's configuration sets, and that order is the author's.
    """
    parities = []
    for level in levels:
        if level.parity not in parities:
            parities.append(level.parity)
    return sorted(levels,
                  key=lambda one: (parities.index(one.parity), one.e_calc))


def header(components=5):
    return FIXED + [name % n for n in range(1, components + 1)
                    for name in ('perc%d', 'conf%d', 'term%d')]


def row_of(level):
    """One line of levels_prepub.txt, as column values.

    The energy formats are levels.txt's, so the two files can be laid
    side by side.  `Ee` and `dE_O-C` are blank together: an unfitted
    level has no observed energy, and levels.txt marks that with the
    `e_c` asterisk beside a copy of `Ec`, which would be worse than
    blank in a file headed for publication.
    """
    out = [level.parity,
           '' if level.unfitted else '%12.6f' % level.e_obs,
           '' if level.unfitted else suffix_text(level.group),
           '%12.4f' % level.e_calc,
           '' if level.unfitted else '%.4f' % (level.e_obs - level.e_calc),
           '%4.1f' % level.j,
           '' if level.g is None else '%6.3f' % level.g]
    for component in level.components:
        out.extend(['%d' % component.percent,
                    component.configuration, component.term])
    return out


def write_prepub(path, levels, components=5):
    with open(path, 'w') as out:
        out.write('\t'.join(header(components)) + '\n')
        for level in levels:
            out.write('\t'.join(row_of(level)).rstrip('\t') + '\n')


def read_prepub(path=PREPUB):
    """Read levels_prepub.txt back into `Level`s, for phase 14.

    The header is the specification, as it is in levels.txt: the fixed
    columns are found by name and the components by the
    `percn`/`confn`/`termn` triples, so a file written with more or
    fewer components reads the same way.  `blocks` stays empty -- the
    per-scheme rows are gone by now, which is the point of the file --
    so what comes back carries the merged composition and nothing that
    would let the merge be second-guessed.
    """
    try:
        with open(path) as stream:
            lines = stream.read().splitlines()
    except IOError:
        raise PurityError('there is no %s here.  make_levels_list_prepub.py '
                          'writes it, from levels.txt and a marked '
                          'coupling_choice.txt.' % path)
    if not lines:
        raise PurityError('%s is empty' % path)
    head = lines[0].split('\t')
    at = {name: i for i, name in enumerate(head)}
    for name in FIXED:
        # `suffix` is not required: a file written before the column
        # existed is still a levels_prepub.txt, and its absence means
        # what a blank cell means -- no level is in a `+x` system.
        if name not in at and name != 'suffix':
            raise PurityError('%s has no %r column; it is not a %s'
                              % (path, name, PREPUB))
    triples, n = [], 1
    while ('perc%d' % n) in at:
        triples.append((at['perc%d' % n], at['conf%d' % n], at['term%d' % n]))
        n += 1

    out = []
    for number, line in enumerate(lines[1:], 2):
        if not line.strip():
            continue
        cell = line.split('\t')

        def get(name):
            i = at.get(name)
            return cell[i].strip() if i is not None and i < len(cell) else ''

        e_calc = _number(get('Ec'), path, number, 'Ec')
        # Blank Ee is what this file says instead of levels.txt's `e_c`
        # asterisk beside a copy of Ec; both mean the same unfitted
        # level, and `row_of` blanks it again on the way out.
        observed = get('Ee')
        level = Level(parity=get('par'),
                      e_obs=(_number(observed, path, number, 'Ee')
                             if observed else e_calc),
                      unfitted=not observed,
                      e_calc=e_calc,
                      j=_number(get('J'), path, number, 'J'),
                      g=(float(get('Lande_g')) if get('Lande_g') else None),
                      conf0='', line=number,
                      # `row_of` wrote it as ' +x'; the space is for
                      # Excel and the '+' is the compilation's, so both
                      # come off and what is left is the letter.
                      group=get('suffix').lstrip('+').lower())
        for percent, conf, term in triples:
            if percent >= len(cell) or not cell[percent].strip():
                break
            level.components.append(Component(
                percent=int(cell[percent]),
                configuration=cell[conf].strip() if conf < len(cell) else '',
                term=cell[term].strip() if term < len(cell) else ''))
        out.append(level)
    if not out:
        raise PurityError('%s carries no levels' % path)
    return out


def _number(text, path, number, name):
    try:
        return float(text)
    except ValueError:
        raise PurityError('%s line %d: %r is not a %s'
                          % (path, number, text, name))
