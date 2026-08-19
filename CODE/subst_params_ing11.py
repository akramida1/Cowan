#!/usr/bin/env python3
"""Copy parameter cards from one ING11 deck into another.

The Python replacement for `subst_params_ing11.pl`.  Its purpose is the
one the Perl's own comment gives: an autoionization run is a separate
deck, built around the continuum configurations, and it should be run
with the parameters that were fitted in the bound-state calculation
rather than with Hartree-Fock ones.  This carries them across.

Usage:

    subst_params_ing11 <source ING11> <target ING11> <set>

`<set>` is the parity set of the **source** deck to take parameters
from, 1 or 2.  The target is rewritten in place, with the deck it had
kept as `<target>.bak`.

A card is matched on its full eighteen-character name, which is the
spectrum and the configuration label together -- the same key RCG uses
to tie a parameter card to a configuration card.  A target card that
the source does not name keeps the parameters it already had, and a
source card the target does not name stops the run: it means the two
decks do not describe the same configurations, and going on would put
one configuration's parameters on another.

Values are carried across as the characters they are punched in and are
never re-encoded, so nothing can lose a digit in transit -- the rule
`toolbox/reorder.py` already works to.

The number of parameters has to agree.  On a single-configuration card
a disagreement is an error, because the two decks then disagree about
what the configuration is; on a CI card the Perl warns and leaves the
target's card alone (`subst_params_ing11.pl:165`), and so does this,
because a CI card's parameter count depends on which pairs the two runs
were asked to couple and the two runs are deliberately different.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files, ing11                                # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


class SubstitutionError(Exception):
    """Two decks that do not describe the same configurations."""


def _named(cards):
    """{18-character name: card}, the key RCG matches a card on."""
    out = {}
    for card in cards:
        out.setdefault(card.name, []).append(card)
    return out


def _carry(targets, sources, kind, warnings):
    """Copy every source card onto the target card of the same name."""
    wanted = _named(sources)
    taken, done = set(), 0
    for card in targets:
        matches = wanted.get(card.name)
        if not matches:
            continue
        source = matches[0]
        taken.add(card.name)
        if source.np != card.np:
            if kind != 'CI':
                raise SubstitutionError(
                    '%s card %s has %d parameters in the source deck and %d '
                    'in the target' % (kind, card.name.strip(), source.np,
                                       card.np))
            warnings.append(
                'CI card %s has %d parameters in the source deck and %d in '
                'the target; it is left as it was.'
                % (card.name.strip(), source.np, card.np))
            continue
        ing11.take_card(card, source)
        done += 1
    missing = sorted(set(wanted) - taken)
    if missing:
        raise SubstitutionError(
            'the target deck has no %s card for %s'
            % (kind, ', '.join(name.strip() for name in missing)))
    return done


def substitute(target, source, kset, warnings):
    """Carry one parity set of `source` onto whatever `target` names it.

    The target's cards are looked at whichever parity set they are in.
    The Perl walked the two files line by line and so needed the target
    to hold a single parity (`subst_params_ing11.pl:6-7`); matching by
    name does not, and a name belongs to one parity anyway.
    """
    if kset > len(source.parameters) or not source.parameters[kset - 1]:
        raise SubstitutionError(
            'the source deck has no parameter card in parity set %d' % kset)
    targets = target.parameters[0] + target.parameters[1]
    target_ci = target.ci[0] + target.ci[1]
    done = _carry(targets, source.parameters[kset - 1], 'configuration',
                  warnings)
    done += _carry(target_ci, source.ci[kset - 1], 'CI', warnings)
    return done


def main(argv):
    if len(argv) != 3:
        raise Usage('three arguments are needed, a source deck, a target '
                    'deck and a parity set')
    source_path, target_path, text = argv
    if text not in ('1', '2'):
        raise Usage('%r is not a parity set of the source deck; use 1 or 2'
                    % text)

    source = ing11.read_ing11(source_path)
    target = ing11.read_ing11(target_path)
    ending = files.line_ending(target_path)

    warnings = []
    done = substitute(target, source, int(text), warnings)
    for warning in warnings:
        sys.stderr.write('subst_params_ing11: %s\n' % warning)
    backup = os.path.basename(target_path) + '.bak'
    files.write_deck(target_path, ing11.write_ing11(target), ending,
                     backup=backup)
    print('%s: %d cards taken from parity set %s of %s; the old deck is %s.'
          % (target_path, done, text, source_path, backup))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('subst_params_ing11: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (SubstitutionError, ing11.Ing11Error, OSError) as problem:
        sys.stderr.write('subst_params_ing11: %s\n' % problem)
        sys.exit(1)
