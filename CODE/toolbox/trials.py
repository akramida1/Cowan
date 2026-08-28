"""The Monte-Carlo machinery behind `calc_TP`.

An uncertainty on a forbidden-line transition probability is not
something Cowan's codes compute.  It is obtained by running RCG many
times over, each run with the fitted Slater parameters and the E2
radial integrals displaced by a random amount drawn from their own
standard deviations, and looking at the spread of the answers.  This
module is the statistics of that: the deviates, the running sums, and
the two Excel-defined shape statistics the author prints to check that
the deviates came out normal.

Everything is accumulated in **ln A**, not in A.  A transition
probability is a squared matrix element and its distribution is far
closer to log-normal than to normal; the spread is therefore reported
as a percentage, `exp(sigma_lnA) - 1`, which is what a log-normal
distribution's relative width is.

THE RANDOM DEVIATES

`calc_TP_M1E2.pl` draws them by hand: two uniforms, a Box-Muller pair
of standard normals, and a scale.  Because Box-Muller yields two
independent normals from one pair of uniforms, the Perl runs the trials
two at a time and requires an even number of them -- the two members of
a pair go to two consecutive trial directories.  `random.gauss` is the
same transformation with the same guarantee, so the pairing has no
purpose here and the trial count may be any number at all.

The generator is seeded explicitly where the parameter file asks for a
seed, which the Perl could not do: a Monte-Carlo run that cannot be
repeated cannot be checked, and these runs take hours.
"""

import math
import random


class Randomizer:
    """Normal deviates, from one seeded stream."""

    def __init__(self, seed=None):
        self.seed = seed
        self.rng = random.Random(seed)

    def normal(self, mean, sigma):
        """One deviate; a zero or negative sigma gives the mean back."""
        if not sigma:
            return mean
        return self.rng.gauss(mean, abs(sigma))


# --------------------------------------------------------------------------
# running sums
# --------------------------------------------------------------------------

class Accumulator:
    """Count, sum and sum of squares of one quantity over the trials."""

    __slots__ = ('n', 'total', 'squares')

    def __init__(self):
        self.n = 0
        self.total = 0.0
        self.squares = 0.0

    def add(self, value):
        self.n += 1
        self.total += value
        self.squares += value * value

    @property
    def mean(self):
        return self.total / self.n if self.n else 0.0

    @property
    def sigma(self):
        """The population standard deviation; never negative from rounding."""
        if not self.n:
            return 0.0
        variance = self.squares / self.n - self.mean * self.mean
        return math.sqrt(variance) if variance > 0 else 0.0


# --------------------------------------------------------------------------
# the two statistics the Perl prints
# --------------------------------------------------------------------------

def stdev(reference, data):
    """The root-mean-square deviation of `data` from `reference`.

    `calc_TP_M1E2.pl:1076-1094`: the population variance about the
    data's own mean, plus the square of the distance from that mean to
    the reference value.  The reference is the unvaried calculation, so
    this is the spread of the trials *and* their bias together.  None
    where there is no data.
    """
    data = list(data)
    if not data:
        return None
    n = len(data)
    average = sum(data) / n
    variance = sum(x * x for x in data) / n - average * average
    if variance < 0:
        variance = 0.0
    offset = reference - average
    return math.sqrt(variance + offset * offset)


def kurtosis(mean, data):
    """Excel's KURT, about a given mean (`calc_TP_M1E2.pl:1048-1073`).

    Printed only as a check on the deviates: a normal distribution has
    kurtosis 0 by this definition, so a number far from 0 says the
    sampling has gone wrong.  None where there is too little data or no
    spread at all.
    """
    data = list(data)
    n = len(data)
    if n < 4:
        return None
    variance = sum(x * x for x in data) / n - mean * mean
    if variance <= 0:
        return None
    total = 0.0
    for x in data:
        d = (x - mean) * (x - mean) / variance
        total += d * d
    return ((total * n * (n + 1) / (n - 1) - 3.0 * (n - 1) * (n - 1))
            / ((n - 2) * (n - 3)))


# --------------------------------------------------------------------------
# the histograms
# --------------------------------------------------------------------------

def rounded_key(x):
    """`calc_TP_M1E2.pl:1035-1044` -- round away from zero, keeping sign.

    Perl's `use integer` truncates towards zero, so the Perl adds 0.5 to
    the magnitude and puts the sign back; `int()` truncates the same way
    and this is the same expression.
    """
    return (-1 if x < 0 else 1) * int(abs(x) + 0.5)


class Histogram:
    """Counts of a deviation in units of one tenth of a sigma."""

    def __init__(self, name):
        self.name = name
        self.bins = {}
        self.values = []

    def add(self, deviation, keep=False):
        """Count one deviation, in sigmas; `keep` also stores the value."""
        self.bins[rounded_key(10.0 * deviation)] = \
            self.bins.get(rounded_key(10.0 * deviation), 0) + 1
        if keep:
            self.values.append(deviation)

    def rows(self):
        """(deviation in sigmas, count), in order."""
        return [(0.1 * k, self.bins[k]) for k in sorted(self.bins)]

    def __bool__(self):
        return bool(self.bins)
