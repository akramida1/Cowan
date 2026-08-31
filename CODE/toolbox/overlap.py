"""Pairing one RCG run's eigenvectors with another's, by overlap.

`levels.identify` ties RCE's fitted levels to RCG's eigenvectors, and it
is right for the unvaried run: it matches on LS composition below the
highest observed level and falls back to position above it, where RCE
did no re-sorting and position is all there is.

That fallback is what fails in a Monte-Carlo trial.  A trial moves the
Slater parameters, the eigenvalues move with them, and two neighbouring
computed levels of the same J change places -- whereupon a positional
match hands the same eigenvector to two different fitted levels and the
run stops.  Measured on Pr III: trial 2 of 5 died with "RCE level 62
(set 2, J=11/2, E=138.961) falls on RCG eigenvector 62, which is already
taken".

The way out is not to identify the trial independently at all.  Both
runs are built from the same deck, so both diagonalize in the same LS
basis, in the same order, with the same basis-state numbering.  A trial
eigenvector is therefore comparable to a base eigenvector directly:

    overlap(i, k) = |<base_i | trial_k>| = |sum_s C_base[s,i] C_trial[s,k]|

one matrix product per (parity set, J) block.  Pair the two runs on
that, largest overlaps first, and compose the result with the base run's
identification.  Every fitted level then keeps the eigenvector it was
identified with, tracked through whatever reordering the trial did.

RCG prints amplitudes to five decimals, so an overlap is good to about
1e-4 -- far finer than the gap between a match (0.99) and a swap.

numpy is a hard requirement here.  `make_cfp`'s optional-numpy
behaviour is unaffected; this module is imported only by `calc_TP`.
"""

try:
    import numpy
except ImportError:                                     # pragma: no cover
    raise ImportError(
        'calc_TP needs numpy to pair the eigenvectors of a trial with '
        'those of the unvaried run.  Install it with:\n'
        '    pip install numpy')


class OverlapError(Exception):
    """Two runs whose eigenvectors cannot be compared."""


def matrix(block, coupling='LS'):
    """One J block's eigenvectors as an array `[basis state, level]`.

    A row of an OUTG11 eigenvector section is one basis state and its
    amplitude in every level of the block, and RCG stops a row at the
    last level it has, so short rows are padded.
    """
    section = block.section(coupling)
    nlevels = block.nlevels
    out = numpy.zeros((len(section.rows), nlevels))
    for at, row in enumerate(section.rows):
        values = row.amplitudes[:nlevels]
        out[at, :len(values)] = values
    return out


def check(base_blocks, trial_blocks):
    """Both runs must have the same blocks and the same basis sizes."""
    bad = []
    for key, block in sorted(base_blocks.items()):
        other = trial_blocks.get(key)
        if other is None:
            bad.append('set %d J=%s is missing from the trial'
                       % (key[0], _j(key[1])))
        elif block.nlevels != other.nlevels:
            bad.append('set %d J=%s has %d levels in the unvaried run and '
                       '%d in the trial'
                       % (key[0], _j(key[1]), block.nlevels, other.nlevels))
    if bad:
        raise OverlapError(
            'the trial does not share a basis with the unvaried run, so '
            'their eigenvectors cannot be compared:\n  '
            + '\n  '.join(bad[:8])
            + ('\n  ... and %d more' % (len(bad) - 8) if len(bad) > 8 else '')
            + '\nThe trial deck is written from the unvaried one, so this '
              'means the trial RCG run did not finish; look at its OUTG11.')


def pair_block(base, trial, coupling='LS'):
    """`([(base level, trial level, overlap)], unpaired)` for one block.

    Levels are 1-based, as RCG numbers them.  The largest overlap is
    taken first and both of its levels struck out, then the next, which
    is the assignment a permutation matrix with a little noise on it
    calls for and is cheap enough at these sizes.
    """
    a, b = matrix(base, coupling), matrix(trial, coupling)
    over = numpy.abs(a.T @ b)                   # [base level, trial level]
    size = over.shape[0]

    order = numpy.argsort(over, axis=None)[::-1]
    taken_base = numpy.zeros(size, dtype=bool)
    taken_trial = numpy.zeros(over.shape[1], dtype=bool)
    pairs = []
    for flat in order:
        i, k = divmod(int(flat), over.shape[1])
        if taken_base[i] or taken_trial[k]:
            continue
        taken_base[i] = taken_trial[k] = True
        pairs.append((i + 1, k + 1, float(over[i, k])))
        if len(pairs) == min(over.shape):
            break
    return pairs, size - len(pairs)


class Report:
    """What the pairing of one whole run looked like."""

    def __init__(self):
        self.overlaps = []
        self.unpaired = 0

    def add(self, pairs, unpaired):
        self.overlaps.extend(value for _, _, value in pairs)
        self.unpaired += unpaired

    @property
    def paired(self):
        return len(self.overlaps)

    @property
    def median(self):
        return float(numpy.median(self.overlaps)) if self.overlaps else 0.0

    def below(self, limit):
        return sum(1 for value in self.overlaps if value < limit)

    def __str__(self):
        text = ('%d eigenvectors paired, median overlap %.3f, %d below 0.7, '
                '%d below 0.3'
                % (self.paired, self.median, self.below(0.7),
                   self.below(0.3)))
        if self.unpaired:
            text += ', %d unpaired' % self.unpaired
        return text


def remap(mapping, base_blocks, trial_blocks, coupling='LS'):
    """Carry the base run's identification over to a trial.

    `mapping` is `{(kset, J, eigenvector): Identified}` for the unvaried
    run, as `calc_TP` builds it.  The return is the same dictionary
    keyed by the trial's eigenvector numbers, and a `Report`.
    """
    check(base_blocks, trial_blocks)
    out = {}
    report = Report()
    for (kset, j), block in sorted(base_blocks.items()):
        pairs, unpaired = pair_block(block, trial_blocks[(kset, j)], coupling)
        report.add(pairs, unpaired)
        for base_level, trial_level, _ in pairs:
            one = mapping.get((kset, j, base_level))
            if one is not None:
                out[(kset, j, trial_level)] = one
    return out, report


def _j(j):
    return '%g' % j if j == int(j) else '%d/2' % round(2 * j)
