"""Ordering a run's configurations the way a published list is ordered.

`conf_list.pl` prints one line per parity, each a comma-separated list of
the run's configurations.  The order is not RCG's and not IN36's: it is
the order such a list is read in, which puts the configuration whose
outermost shell is lowest first, and among those the one whose inner
shells are fullest.

The sort key below is the Perl's, kept exactly, because the order is a
presentation choice the author made and not something derivable:

    (1000 + n of the outermost occupied subshell)
    then, for each subshell up to that one, (n x 10000, l, 99 - occupation)

compared as **text**, which is what `conf_list.pl:37` does.  Padding `n`
to four digits and inverting the occupation are what make a text compare
give the intended numeric order, and `99 - occupation` is what puts the
fuller subshell first.

What the subshells are comes from `toolbox/in36.py` -- RCN's own ANALYZ
and SETCFG -- rather than from the Perl's fixed-offset regex.  Only the
subshells the *card* names are used, not the noble-gas expansion, so the
printed list is the one the author punched.
"""

from .in36 import LBCD

# The Perl's own l table, `conf_list.pl:6-7`.  It stops at l = 11 where
# RCN's LBCD runs to l = 20, and nothing in these decks goes past 'o'.
_MAX_L = 11


class ConfListError(Exception):
    """A configuration that cannot be put in a printed list."""


def card_shells(configuration):
    """The subshells the IN36 card names, in card order.

    `conf_list.pl` reads at most eight of them (`:41`, `:117`), because
    Cowan's card holds eight, and so does this.
    """
    return configuration.card_shells[:8]


def drop_leading_f14(configurations):
    """Whether every configuration opens on a closed f shell.

    `conf_list.pl:94-105` drops the first subshell of every
    configuration when all of them begin with f14.  It is the one closed
    shell these decks routinely carry and never want printed; any other
    shared closed shell is left in, because dropping it would change the
    electron count of the printed labels and the plan forbids that
    except at the user's explicit request.
    """
    shells = [card_shells(c) for c in configurations]
    if not shells or not all(s for s in shells):
        return False
    return all(s[0].l == 3 and s[0].occ == 14 for s in shells)


def _shells_of(configuration, drop_first):
    shells = card_shells(configuration)
    return shells[1:] if drop_first else shells


def _last_occupied(shells):
    """The index of the outermost subshell that holds electrons.

    `get_last_shell` (`conf_list.pl:112-126`) scans from the end and
    returns 0 when the configuration is empty, which is the same as
    stopping at the first subshell; that is kept.
    """
    for i in range(len(shells) - 1, -1, -1):
        if shells[i].occ > 0:
            return i
    return 0


def sort_key(configuration, drop_first=False):
    """The Perl's text sort key, `conf_list.pl:16-28`."""
    shells = _shells_of(configuration, drop_first)
    if not shells:
        return ''
    last = _last_occupied(shells)
    key = ''
    n_max = 0
    for shell in shells[:last + 1]:
        if shell.l > _MAX_L:
            raise ConfListError(
                'configuration %s has an %s subshell, which is past the '
                'l table conf_list uses' % (configuration.label,
                                            LBCD[shell.l]))
        key += '%d %d %d' % (shell.n * 10000, shell.l, 99 - int(shell.occ))
        n_max = max(n_max, shell.n)
    return '%d %s' % (1000 + n_max, key)


def parity_of(configuration, drop_first=False):
    """0 even, 1 odd, from the subshells that are printed.

    Taken from the printed subshells rather than from
    `Configuration.parity` so that dropping a closed f shell cannot
    change it -- a closed shell is even, so the two always agree, and
    this way the label and the parity it is filed under are built from
    the same list.
    """
    odd = sum(int(s.occ) for s in _shells_of(configuration, drop_first)
              if s.l % 2)
    return odd % 2


def label(configuration, drop_first=False):
    """One configuration in ordinary notation: `1s2.2s`, `4f.5p6.5d2`.

    Empty subshells are left out and an occupation of one is not
    written, which is `conf_list.pl:41-52`.
    """
    out = []
    for shell in _shells_of(configuration, drop_first):
        if shell.occ:
            occ = '' if shell.occ == 1 else '%d' % int(shell.occ)
            out.append('%s%s' % (shell.nl, occ))
    return '.'.join(out)


def parity_line(configurations, parity):
    """The comma-separated list of one parity, in the printed order."""
    drop_first = drop_leading_f14(configurations)
    wanted = [c for c in configurations
              if parity_of(c, drop_first) == parity]
    wanted.sort(key=lambda c: sort_key(c, drop_first))
    return ', '.join(label(c, drop_first) for c in wanted)
