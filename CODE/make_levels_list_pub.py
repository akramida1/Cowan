#!/usr/bin/env python3
"""Write levels_pub.txt: a designation no other level of its (parity, J) has.

The last of the three tools that turn a finished fit into a published
level list:

    levels.txt --12--> coupling_choice.txt --13--> levels_prepub.txt
                                                   --14--> levels_pub.txt

`make_levels_list_prepub.py` wrote each configuration's components in
the coupling scheme chosen for it.  This names the levels.  A published
or database level needs a name no other level of the same parity and J
carries: the leading component supplies it where the level is
reasonably pure, and where it is not the name comes from further down
the composition, then from letter prefixes, then from the energy.  See
`toolbox/publabel.py` for why that order and what each step costs.

Usage:

    make_levels_list_pub [--prepub levels_prepub.txt] [-o levels_pub.txt]
                         [--nlev 555,481]        levels per configuration set
                         [--n_comp 3]            composition groups, 1 to 5
                         [--min-name-percent 20] smallest component to name by
                         [--letters-max 50]      % of levels allowed letters
                         [--leading-min 50]      % allowed under a 50% leading
                         [--max-outside 5]       % named from beyond n_comp

`--nlev` is one count per configuration set in IN36 order -- in practice
one per parity, RCE taking two sets and no more.  Its default is the
position of the highest observed level in each set; a set with no
observed level makes it mandatory.

The energies are written in cm-1.  Everything above this file carries
RCE's own 10^3 cm-1 and is left in it; this is the file that gets
published, so it is the one converted.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import prepub, publabel                           # noqa: E402
from toolbox.purity import PurityError                         # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


class Options:

    def __init__(self, argv):
        self.prepub = prepub.PREPUB
        self.out = publabel.PUB
        self.nlev = None
        self.n_comp = 3
        self.min_name_percent = 20.0
        self.letters_max = 50.0
        self.leading_min = 50.0
        self.max_outside = 5.0

        rest = list(argv)
        positional = []
        while rest:
            arg = rest.pop(0)
            if arg == '-o':
                if not rest:
                    raise Usage('-o needs a file name')
                self.out = rest.pop(0)
                continue
            if not arg.startswith('--'):
                positional.append(arg)
                continue
            name, _, value = arg[2:].partition('=')
            name = name.replace('-', '_')
            if not value:
                if not rest:
                    raise Usage('--%s needs a value' % name.replace('_', '-'))
                value = rest.pop(0)
            if name == 'prepub':
                self.prepub = value
            elif name in ('o', 'out'):
                self.out = value
            elif name == 'nlev':
                self.nlev = _counts(value)
            elif name == 'n_comp':
                self.n_comp = _number(int, value, '--n_comp')
            elif name == 'min_name_percent':
                self.min_name_percent = _number(float, value,
                                                '--min-name-percent')
            elif name == 'letters_max':
                self.letters_max = _number(float, value, '--letters-max')
            elif name == 'leading_min':
                self.leading_min = _number(float, value, '--leading-min')
            elif name == 'max_outside':
                self.max_outside = _number(float, value, '--max-outside')
            else:
                raise Usage('--%s is not an option' % name.replace('_', '-'))

        if positional:
            self.prepub = positional[0]
        if len(positional) > 1:
            raise Usage('too many file names: %s' % ', '.join(positional))
        if not 1 <= self.n_comp <= 5:
            raise Usage('--n_comp is between 1 and 5, not %d' % self.n_comp)


def _number(kind, text, name):
    try:
        return kind(text)
    except ValueError:
        raise Usage('%s wants a number, not %r' % (name, text))


def _counts(text):
    """`--nlev 555,481`: one count per configuration set, in IN36 order."""
    out = []
    for piece in text.replace(',', ' ').split():
        value = _number(int, piece, '--nlev')
        if value < 1:
            raise Usage('--nlev counts at least one level per set, not %r'
                        % piece)
        out.append(value)
    if not out:
        raise Usage('--nlev needs a count for each configuration set')
    return out


def main(argv):
    options = Options(argv)
    levels = prepub.read_prepub(options.prepub)
    published, counts = publabel.cut(levels, options.nlev)
    print('%s: %d levels; publishing %d (%s)%s'
          % (options.prepub, len(levels), len(published),
             ', '.join('%s %d' % (parity, counts[parity])
                       for parity in publabel.sets(levels)),
             '' if options.nlev else
             ', to the highest observed level of each set'))

    chosen = publabel.choose(published,
                             min_percent=options.min_name_percent)
    names, lettered, nameless = publabel.name_all(chosen)

    why, share = publabel.fallback(names, lettered,
                                   letters_max=options.letters_max,
                                   leading_min=options.leading_min)
    if why:
        indexed = publabel.by_energy(names)
        print('   %s, so a component cannot name these levels at all and '
              'every one of them is named by %s.  The composition beside it '
              'is what carries the physics.'
              % ('%.0f%% of the published levels need a letter, over the '
                 '--letters-max of %g%%' % (share, options.letters_max)
                 if why == 'letters' else
                 '%.0f%% of the published levels have a leading component '
                 'under 50%%, over the --leading-min of %g%%'
                 % (share, options.leading_min),
                 'a sequential index within its parity, two integer energies '
                 'of one parity being equal' if indexed
                 else 'the integer part of its energy in cm-1'))
    else:
        _report(options, names, lettered, nameless)

    clash = publabel.duplicates(names)
    if clash:
        first, second = clash[0]
        raise publabel.PubError(
            '%d designations are not unique within their (parity, J); the '
            'first is %s %s at Ec=%.3f and Ec=%.3f, both %s J=%.1f.  That is '
            'a defect in the search, not in the data.'
            % (len(clash), first.label[0], first.label[1],
               first.level.e_calc * publabel.CM1,
               second.level.e_calc * publabel.CM1,
               first.level.parity, first.level.j))

    publabel.write_pub(options.out, names, n_comp=options.n_comp)
    print('wrote %s (%d levels, %d composition group%s)'
          % (options.out, len(names), options.n_comp,
             '' if options.n_comp == 1 else 's'))
    return 0


def _report(options, names, lettered, nameless):
    """The three things the author asked to be told about by name."""
    total = len(names) or 1
    if lettered:
        print('   %d levels (%.0f%%) share a name with another of their '
              'parity and J and take a letter prefix: %s%s'
              % (len(lettered), 100.0 * len(lettered) / total,
                 ', '.join('%s %s (%s Ec=%.1f J=%.1f)'
                           % (one.label[1], one.conf, one.level.parity,
                              one.level.e_calc * publabel.CM1, one.level.j)
                           for one in lettered[:5]),
                 ' ...' if len(lettered) > 5 else ''))
    if nameless:
        print('   %d levels have no component left to be named after and are '
              'named by their energy: %s'
              % (len(nameless),
                 ', '.join('%s J=%.1f' % (one.level.parity, one.level.j)
                           for one in nameless[:5])))

    # The author's ruling: a name may come from outside the printed
    # groups, and that is accepted -- but a file where it happens often
    # is a file whose compositions do not explain their own names, so
    # the levels are listed once it passes --max-outside.
    beyond = publabel.outside(names, options.n_comp)
    if beyond:
        share = 100.0 * len(beyond) / total
        print('   %d levels (%.0f%%) are named after a component past the '
              '%d printed%s'
              % (len(beyond), share, options.n_comp,
                 ':' if share > options.max_outside
                 else ', within the --max-outside of %g%%'
                      % options.max_outside))
        if share > options.max_outside:
            for one in beyond:
                print('      %s Ec=%12.3f J=%4.1f  %d%% %s %s (component %d)'
                      % (one.level.parity, one.level.e_calc * publabel.CM1,
                         one.level.j,
                         one.percent, one.conf, one.label[1], one.index + 1))


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('make_levels_list_pub: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except PurityError as problem:
        sys.stderr.write('make_levels_list_pub: %s\n' % problem)
        sys.exit(1)
