#!/usr/bin/env python3
"""Print the configurations of IN36 as one sorted list per parity.

The Python replacement for `conf_list.pl`.  The output is two lines,
even parity then odd, each a comma-separated list of configurations in
ordinary spectroscopic notation:

    1s2.2s, 1s2.2p2, 1s.2s2, ...
    1s2.2p, 1s.2s.2p

It is a convenience for pasting a configuration list into a paper or a
database entry, so the order is the reading order of such a list --
outermost shell last, lowest configuration first -- and not RCG's.

Usage:

    conf_list [folder]

IN36 is read from `folder`, or from the working directory.

What this does that the Perl did not: it parses IN36 through
`toolbox/in36.py`, which is RCN's own ANALYZ/SETCFG, instead of the one
fixed-offset regex at `conf_list.pl:68`.  That regex reads a subshell
only where RCN's free format happens to put it, and it is one of the
three independent IN36 parsers the plan set out to remove.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import conflist, in36                              # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


def main(argv):
    if len(argv) > 1:
        raise Usage('at most one argument is allowed, the folder to read')
    folder = argv[0] if argv else '.'
    configurations = in36.read_in36(folder=folder)
    for parity in (0, 1):
        print(conflist.parity_line(configurations, parity))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('conf_list: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (in36.In36Error, OSError) as problem:
        sys.stderr.write('conf_list: %s\n' % problem)
        sys.exit(1)
