"""Autoionization rates, as `read_aa` tabulates them.

One file per parity set, tab separated, written by `read_aa.pl` (phase 9
turns it into `read_aa.py`) out of the Auger widths RCG prints:

    J   Energy   1s2:0   1s2s:120.1   ...   AA_tot

The header names one column per autoionization channel, each with the
energy of the core state it leaves behind, in the run's own energy unit
and relative to this ion's ground state.  The rows that follow are the
levels of one parity set, in ascending J and, within a J, in the order
of RCG's eigenvectors -- so a row is tied to a level by position, and
the J it carries is checked against the level's rather than trusted.

The first channel carries the dielectronic recombination rate: it is
the one to the ground state of the ionized atom, which is what a DR
spectrum's population follows (`conv_out.pl:ReadAI`).
"""

from dataclasses import dataclass


class AiError(Exception):
    """An autoionization file that does not match the run."""


@dataclass(frozen=True)
class Channel:
    """One autoionization channel: the core it leaves, and that core's energy."""

    name: str
    energy: float

    @property
    def text(self):
        """The column heading, which is how `conv_out` labels the channel."""
        return '%s:%s' % (self.name, _plain(self.energy))


@dataclass
class Rates:
    """One level's rates, one per channel, plus the total."""

    j: float
    energy: float
    rates: dict                 # channel name -> rate, s-1
    total: float

    def dr(self, channels):
        """The rate to the ground state of the ionized atom.

        That is the first channel, the one `read_aa` is told to put
        first, and a dielectronic-recombination spectrum's populations
        follow it.
        """
        return self.rates.get(channels[0].name, 0.0) if channels else 0.0


def _plain(value):
    return '%g' % value


def read_ai(path):
    """Read one `read_aa` table."""
    with open(path, 'r', errors='replace') as fh:
        text = fh.read()
    rows = [line.rstrip('\r').split('\t')
            for line in text.split('\n') if line.strip()]
    if not rows:
        raise AiError('%s is empty' % path)
    heading = rows[0]
    if len(heading) < 4 or heading[0].strip() != 'J':
        raise AiError(
            '%s does not start with the heading read_aa writes '
            '(J, Energy, one column per channel, AA_tot)' % path)
    channels = []
    for column in heading[2:-1]:
        name, _, energy = column.strip().partition(':')
        try:
            channels.append(Channel(name=name, energy=float(energy)))
        except ValueError:
            raise AiError('%r is not a channel heading of the form core:energy'
                          % column)

    out = []
    for at, row in enumerate(rows[1:], 2):
        if len(row) != len(heading):
            raise AiError('%s line %d has %d columns and its heading %d'
                          % (path, at, len(row), len(heading)))
        rates = {}
        for channel, cell in zip(channels, row[2:-1]):
            value = float(cell or 0)
            if value > 0:
                rates[channel.name] = value
        out.append(Rates(j=float(row[0]), energy=float(row[1]), rates=rates,
                         total=float(row[-1] or 0)))
    return channels, out


def pair(rows, identified, blocks, path):
    """Tie one file's rows to the levels of a parity set, by position.

    The order is ascending J and then RCG eigenvector number, which is
    the order `read_aa` writes; the J on the row must agree with the
    level's or the two files are not of the same run.
    """
    order = sorted(identified, key=lambda one: (one.j, one.rcg))
    if len(order) != len(rows):
        raise AiError('%s has %d levels and the parity set %d'
                      % (path, len(rows), len(order)))
    out = {}
    for one, row in zip(order, rows):
        if abs(one.j - row.j) > 1e-6:
            raise AiError(
                'J disagrees between %s and OUTG11: the file says %s where '
                'the run has %s, at E = %s' % (path, row.j, one.j, row.energy))
        out[(one.kset, one.j, one.rcg)] = row
    return out
