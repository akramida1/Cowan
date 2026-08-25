#!/usr/bin/env python3
"""Write levels.txt and params.txt from a finished RCG/RCE run.

The Python replacement for `printout.pl`.  What it does that the Perl
could not:

* **every coupling RCG can print**, not LS and Cowan's jj alone, named
  on the command line and, where the run that has one is a second
  folder, taken from there (`--from LSJK=JK`).  That merge was hand
  work; see `toolbox/levels.py` for the two ways it can be done and how
  the right one is chosen.
* **labels built from `genealogy.txt`**, so the middle of a coupling
  chain is read rather than guessed, and **n taken from IN36 through the
  ING11 slots**, so a reordered deck is labelled correctly.
* **one electron count in every label.**  Only `--drop-full-shells` can
  break that, and only where the named subshell is closed.

Usage:

    printout <levels file> [<parameters file>]
             [--couplings LS,JJ]        what to print, in this order
             [--from LSJK=JK]           take a coupling from a sub-folder
             [--from LSr=LS@reordered]  ... under a name of your own
             [--dominant-config]        add the conf0 column
             [--drop-full-shells 5p,4d] drop these where they are closed
             [--match eigenvector|energy-order]
             [--components 5] [--min-percent 1.5]
             [--no-rce]                 RCG only, no fitted energies

IN36 is read from the folder above where the run's own folder has none,
which is what the `JK/` arrangement looks like from the inside.  --no-rce
is assumed where there is no RCEOUT to read, so a folder RCG has run in
and RCE has not needs no option at all.

Defaults are the Perl's: levels.txt, params.txt, five components, 1.5%.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import case as case_mod                            # noqa: E402
from toolbox import asd, designate, files, in36, levels, rce   # noqa: E402


def header(options):
    """The columns of levels.txt.

    `conf0` is off by default, and with it off the file is exactly what
    it has always been: nothing downstream of this script but the
    phase 12-14 pipeline wants the column, and levels.txt is read by
    hand and by spreadsheet as well as by code.
    """
    out = ['cpl', 'par', 'Ee', 'e_c', 'Ec', 'Ee-Ec', 'J', 'Lande_g']
    if options.dominant_config:
        out.append('conf0')
    # As many triples as `--components` asks for, not the Perl's fixed
    # five: the header is what tells a reader -- and `purity.py` --
    # where the components are, so five names over ten written triples
    # would silently hide the last five.  With the default it is the
    # same header it has always been.
    return out + [c % n for n in range(1, options.components + 1)
                  for c in ('%%%d', 'conf%d', 'term%d')]


PARAM_HEADER = ['parity', 'Configuration', '', 'parameter', 'LSF', 'Uncert.',
                'Group', 'HF', 'LSF/HF']

#: Which HR scale factor divides which kind of parameter back out of the
#: value RCN2 punched: Fk(i,i), zeta, Fk(i,j), Gk(i,j), and the fifth for
#: the configuration-interaction integrals (`rcg11kd.f:10837`, `:10905`).
_SCALE_OF = {'F': 0, 'ZETA': 1, 'f': 2, 'G': 3, 'R': 4}


class Usage(Exception):
    """The command line will not do what was asked."""


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------

class Block:
    """One block of the level list: a column name, what it is, where from.

    `name` is what goes in the `cpl` column and is normally the coupling
    itself.  It is separate from `coupling` only so that one coupling
    can be printed twice from two folders -- `--from LSr=LS@reordered`
    -- which a single levels.txt otherwise cannot hold, there being one
    `cpl` value per block.  The alias is the author's to choose: nothing
    downstream may invent a name for a block.
    """

    def __init__(self, name, coupling, folder=None):
        self.name = name
        self.coupling = coupling
        self.folder = folder


class Options:

    def __init__(self, argv):
        self.levels = 'levels.txt'
        self.params = 'params.txt'
        self.couplings = ['LS']
        self.sources = {}               # block name -> (coupling, folder)
        self.blocks = []                # [Block], in print order
        self.drop = []                  # (n, l)
        self.match = None
        self.components = 5
        self.min_percent = 1.5
        self.no_rce = False
        self.dominant_config = False

        positional = []
        rest = list(argv)
        while rest:
            arg = rest.pop(0)
            if not arg.startswith('--'):
                if arg.lower().rstrip('=') in ('norce', 'no_rce'):
                    self.no_rce = True
                else:
                    positional.append(arg)
                continue
            name, _, value = arg[2:].partition('=')
            name = name.replace('-', '_')
            if name == 'no_rce':
                self.no_rce = True
                continue
            if name == 'dominant_config':
                self.dominant_config = True
                continue
            if not value:
                if not rest:
                    raise Usage('--%s needs a value' % name.replace('_', '-'))
                value = rest.pop(0)
            if name == 'couplings':
                # Not validated here: a name may be an alias that a
                # --from still to be read will define, and the options
                # may come in either order.
                self.couplings = [c.strip() for c in value.split(',')
                                  if c.strip()]
            elif name == 'from':
                block, _, where = value.partition('=')
                if not where:
                    raise Usage('--from wants COUPLING=FOLDER, not %r' % value)
                block = block.strip()
                coupling, _, folder = where.partition('@')
                if folder:
                    # The aliased form: a name of the author's, the
                    # coupling it holds, and the folder it comes from.
                    self.sources[block] = (_coupling(coupling), folder)
                else:
                    self.sources[block] = (_coupling(block), coupling)
            elif name == 'drop_full_shells':
                self.drop = [_nl(s) for s in value.split(',') if s.strip()]
            elif name == 'match':
                if value not in (levels.EIGENVECTOR, levels.ENERGY_ORDER):
                    raise Usage('--match is %s or %s, not %r'
                                % (levels.EIGENVECTOR, levels.ENERGY_ORDER,
                                   value))
                self.match = value
            elif name == 'components':
                self.components = int(value)
            elif name == 'min_percent':
                self.min_percent = float(value)
            else:
                raise Usage('--%s is not an option' % name.replace('_', '-'))

        if positional:
            self.levels = positional[0]
        if len(positional) > 1:
            self.params = positional[1]
        if len(positional) > 2:
            raise Usage('too many file names: %s' % ', '.join(positional))
        unknown = [c for c in self.sources if c not in self.couplings]
        if unknown:
            raise Usage('--from names %s, which --couplings does not ask for'
                        % ', '.join(unknown))
        seen = set()
        for at, name in enumerate(self.couplings):
            if name in self.sources:
                coupling, folder = self.sources[name]
                block = Block(name, coupling, folder)
            else:
                # No --from, so the name has to be a coupling of the
                # primary folder and is validated as one.
                self.couplings[at] = _coupling(name)
                block = Block(self.couplings[at], self.couplings[at])
            if block.name in seen:
                raise Usage('--couplings asks for %s twice; two blocks of one '
                            'coupling need a name each (--from %sr=%s@folder)'
                            % (block.name, block.name, block.name))
            seen.add(block.name)
            self.blocks.append(block)


def _coupling(name):
    name = name.strip()
    for known in designate.TABLE_OF:
        if known.lower() == name.lower():
            return known
    raise Usage('%r is not one of RCG\'s couplings (%s)'
                % (name, ', '.join(sorted(designate.TABLE_OF))))


def _nl(text):
    """`5p` -> (5, 1)."""
    text = text.strip()
    for i, ch in enumerate(text):
        if ch.isalpha():
            try:
                return (int(text[:i]), in36.LBCD.index(ch.lower()))
            except ValueError:
                break
    raise Usage('%r is not a subshell like 5p' % text)


# --------------------------------------------------------------------------
# the level list
# --------------------------------------------------------------------------

def _percent(amplitude):
    """RCG prints amplitudes; a percentage is the square, rounded."""
    return int(amplitude * amplitude * 100 + 0.5)


def _components(case, coupling, kset, j, eigenvector, options):
    """The largest components of one eigenvector, as report columns."""
    weights = case.amplitudes(coupling, kset, j, eigenvector)
    out = []
    for rank, (amplitude, state) in enumerate(weights[:options.components], 1):
        percent = _percent(amplitude)
        if percent < options.min_percent and rank > 1:
            break
        found = case.designation(coupling, kset, j, state)
        if not found.canonical:
            _NONCANONICAL.add((kset, j, state))
        out.extend(['%d' % percent, found.configuration, found.term])
    return out


#: Basis states printed in the deck's summation order because an
#: intermediate term made that order observable; see `designate.py`.
_NONCANONICAL = set()


def _rows(primary, block, source, options):
    """One coupling's block of the level list, both parities, in RCG
    eigenvalue order, as lists of column values.

    The energies and the identification are always the primary run's;
    only the eigenvector the components are read from moves to the
    second folder, through the twin map the merge produced.  A level
    with no twin keeps its energies and loses its components, which is
    reported rather than silently dropped.
    """
    case = source or primary
    coupling = block.coupling
    twin = getattr(source, 'twin', None) if source else None
    out = []
    lost = 0
    for kset in sorted({k for k, _ in primary.blocks}):
        parity = 'e' if primary.card(kset, 1).parity == 0 else 'o'
        found = [one for one in primary.identified if one.kset == kset]
        if not found:
            found = levels.unfitted(primary.blocks, kset)
        for one in sorted(found,
                          key=lambda x: levels.order_key(primary.blocks, x)):
            level = one.level
            observed = '' if level.observed else '*'
            e_obs = level.e_calc if not level.observed else level.e_obs
            difference = ('' if not level.observed
                          else '%.4f' % (e_obs - level.e_calc))
            row = [block.name, parity,
                   '%12.6f' % e_obs, observed, '%12.4f' % level.e_calc,
                   difference, '%4.1f' % one.j,
                   '' if level.g is None else '%6.3f' % level.g]
            if options.dominant_config:
                row.append(_readable(primary, kset, level.conf0))
            eigenvector = one.rcg
            if twin is not None:
                eigenvector = twin.get((kset, one.j, one.rcg))
            if eigenvector is None:
                lost += 1
            else:
                row.extend(_components(case, coupling, kset, one.j,
                                       eigenvector, options))
            out.append(row)
    if lost:
        print('   %d %s levels have no twin and carry no components'
              % (lost, block.name))
    return out


# --------------------------------------------------------------------------
# the parameter list
# --------------------------------------------------------------------------

def _scale_index(named):
    """Which of the four HR factors applies, or None where none does.

    They are punched in the order Fk(i,i), zeta, Fk(i,j), Gk(i,j), so an
    F integral takes the first or the third according to whether its two
    subshells are the same one (`rcg11kd.f:10837`).  Eav and the
    effective parameters are not scaled at all.
    """
    if named.kind == 'ZETA':
        return _SCALE_OF['ZETA']
    if named.kind == 'G':
        return _SCALE_OF['G']
    if named.kind == 'F':
        return _SCALE_OF['F' if len(set(named.shells)) == 1 else 'f']
    return None


def _param_rows(case, oute, hf, options):
    """params.txt, one row per fitted parameter.

    OUTE and RCEINP.HF list the same parameters in the same order as
    RCEOUT -- and so does OUTGINE, which stands in for both where RCE
    has not run -- and they are paired by position rather than by name: the
    name RCE prints is a slot number, so `ZETA 3` occurs once per
    configuration that has a slot 3, and a lookup by name silently
    collapses them (`conv_cowan.pl` keys them by index too).

    The fitted value comes from RCEOUT, which prints a digit more than
    either of the others; the least-squares flag and the standard
    deviation come from OUTE, which is the only file that has them.
    """
    out = []
    for kset in sorted({k for k, _ in case.blocks}):
        parity = 'e' if case.card(kset, 1).parity == 0 else 'o'
        named = case.parameters(kset)
        fitted = case.parameter_source(kset)
        sigmas = list(oute.get(kset, ()))
        starting = list(hf.get(kset, ()))
        if sigmas and len(sigmas) != len(fitted):
            raise case_mod.CaseError(
                'OUTE lists %d parameters for parity set %d and RCEOUT %d'
                % (len(sigmas), kset, len(fitted)))
        groups = {}
        for one in sigmas or fitted:
            groups[abs(one.flag)] = groups.get(abs(one.flag), 0) + 1

        for at, (one, parameter) in enumerate(zip(named, fitted)):
            sigma = sigmas[at] if at < len(sigmas) else None
            start = starting[at] if at < len(starting) else None
            flag = sigma.flag if sigma is not None else parameter.flag
            hf_value = None
            if start is not None:
                hf_value = start.value / _factor(case, kset, one)
            uncertainty = ''
            if sigma is not None:
                uncertainty = ('fixed' if abs(flag) >= 99
                               else '%.3f' % (abs(sigma.sigma) * 1000))
            group = ''
            if 1 <= abs(flag) < 99 and groups.get(abs(flag), 0) > 1:
                group = '%d' % abs(flag)
            ratio = ('' if not hf_value
                     else '%.8f' % (parameter.value / hf_value))
            row = [parity]
            row.extend(_configuration_columns(case, kset, one))
            row.extend([one.text, '%.3f' % (parameter.value * 1000),
                        uncertainty, group,
                        '' if hf_value is None else '%.3f' % (hf_value * 1000),
                        ratio])
            out.append(row)
    return out


def _factor(case, kset, named):
    """The HR factor a Hartree-Fock value has to be divided by.

    RCN2 punches the parameters into ING11 already multiplied by the
    scale factors on the card (`IFACT`, columns 73-80), and RCE writes
    its starting values out the same way, so a comparison of the fitted
    value with the *unscaled* Hartree-Fock one has to take them off
    again -- which is what the LSF/HF column means.
    """
    if named.kind == 'R':
        return case.ci_factor(kset, named.configurations)
    index = _scale_index(named)
    if index is None:
        return 1.0
    number = 1
    for at, card in enumerate(case.deck.set(kset), 1):
        if card.label in named.configurations:
            number = at
            break
    four = case.scale_factors(kset, number)
    return four[index] if index < 4 else four[3]


def _configuration_columns(case, kset, named):
    """The one or two configuration columns of a parameter row."""
    labels = list(named.configurations) or ['']
    if len(labels) == 1:
        labels.append('')
    return [_readable(case, kset, labels[0]), _readable(case, kset, labels[1])]


def _readable(case, kset, label):
    """`f26p` -> `4f2.5p6.6p`: the configuration as a reader knows it."""
    if not label:
        return ''
    for card in case.deck.set(kset):
        if card.label == label:
            shells = case.shells[(card.kset, card.index)]
            tokens = [asd.Subshell(n=s.n, l=s.l, occ=s.occ)
                      for s in shells if s.occ
                      and not ((s.n, asd._l_value(s.l)) in case.drop
                               and not s.is_open)]
            tokens.sort(key=lambda t: (t.n, t.lval))
            return asd.configuration(tokens)
    return label


# --------------------------------------------------------------------------
# putting it out
# --------------------------------------------------------------------------

def _write(path, header, rows):
    with open(path, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for row in rows:
            out.write('\t'.join(row).rstrip('\t') + '\n')


def main(argv):
    options = Options(argv)
    wanted = []
    for block in options.blocks:
        if block.folder is None and block.coupling not in wanted:
            wanted.append(block.coupling)
    # A folder RCG has run in but RCE has not is a whole case of its own
    # -- the second coupling of a fit, say -- so it is reported on rather
    # than refused.  --no-rce is what that report is, and asking for it
    # explicitly is only a way of skipping OUTE and RCEINP.HF where RCE
    # has in fact run.
    if not options.no_rce and files.resolve(rce.RCEOUT, '.',
                                            required=False) is None:
        options.no_rce = True
        print('no %s here, so reporting the calculated levels alone, as '
              '--no-rce does' % rce.RCEOUT)
    primary = case_mod.Case.read('.', couplings=wanted or ['LS'],
                                 with_rce=not options.no_rce,
                                 drop=options.drop)
    # With --no-rce there is nothing identified against a fit, so the
    # count is RCG's own eigenvalues -- which is what the blocks below
    # will print, and what the reader is counting.
    found = len(primary.identified) or sum(
        len(block.eigenvalues) for block in primary.blocks.values())
    print('%s: %d J blocks, %d levels'
          % (', '.join(wanted) or 'LS', len(primary.blocks), found))

    sources = {}
    for block in options.blocks:
        if block.folder is None:
            continue
        other = case_mod.Case.read(block.folder,
                                   couplings=[block.coupling, 'LS'],
                                   parent='.', with_rce=not options.no_rce,
                                   drop=options.drop)
        mode = options.match or levels.match_mode(primary.deck.configurations,
                                                  other.deck.configurations)
        merged = _merge(primary, other, mode)
        sources[block.name] = _relabel(other, merged)
        print('%s (%s) from %s, matched by %s'
              % (block.name, block.coupling, block.folder, mode))
        if mode == levels.ENERGY_ORDER:
            _report_energy_order(merged)

    rows = []
    for block in options.blocks:
        rows.extend(_rows(primary, block, sources.get(block.name), options))
    _write(options.levels, header(options), rows)
    print('wrote %s (%d levels)' % (options.levels, len(rows)))
    if _NONCANONICAL:
        print('   %d basis states are printed in the deck\'s summation '
              'order, not in n,l order: an intermediate term makes that '
              'order visible, so reordering the subshells would change '
              'what the label says' % len(_NONCANONICAL))
    for coupling, kset, j in sorted(
            set(primary.unlabelled)
            | {u for s in sources.values() for u in s.unlabelled}):
        print('   %s J=%s of set %d has one basis state, so RCG printed no '
              '%s section; it is named in LS' % (coupling, ('%.1f' % j),
                                                 kset, coupling))

    if not options.no_rce:
        oute = rce.read_oute(folder='.')
        hf = rce.read_hf(folder='.')
        parameters = _param_rows(primary, oute, hf, options)
        _write(options.params, PARAM_HEADER, parameters)
        print('wrote %s (%d parameters)' % (options.params, len(parameters)))
    else:
        _write_params_without_rce(primary, options)
    return 0


def _write_params_without_rce(primary, options):
    """params.txt from the deck alone, with the HF column reconstructed.

    Nothing here has been fitted, so the LSF column is the value the run
    itself used -- what ING11 holds and OUTGINE punched.  The HF column
    is that value divided back by the HR scale factor of its card, which
    is what RCEINP.HF would have carried had RCE been run on an
    unfitted deck; `_factor` is the same division the fitted path makes,
    so the LSF/HF column comes out as the scale factor.

    Two things about it are not the real Hartree-Fock value and are said
    so at the end of the run rather than left to be discovered.  Eav is
    not scaled at all, so its HF column is the deck's own number
    repeated.  And where the deck has been through `s11`, the values are
    a fit's and the factors are still the original ones, so dividing one
    by the other recovers nothing: only RCEINP.HF, written by an RCE run
    over the unfitted deck, has those numbers.
    """
    if not primary.starting:
        print('no OUTGINE here, so no %s: it is the only file in a folder '
              'RCE has not run in that carries the parameter names'
              % options.params)
        return
    parameters = _param_rows(primary, {}, primary.starting, options)
    _write(options.params, PARAM_HEADER, parameters)
    print('wrote %s (%d parameters)' % (options.params, len(parameters)))
    print("   --no-rce, so the values are the deck's own and the HF column "
          "is each one divided by the scale factor of its card; Eav has no "
          "scale factor and is repeated as it stands.  RCEINP.HF, from an "
          "RCE run over the unfitted deck, is what makes every HF value "
          "right")


def _merge(primary, other, mode):
    if mode == levels.EIGENVECTOR:
        return levels.merge_by_eigenvector(primary.identified,
                                           other.identified,
                                           primary.blocks, other.blocks)
    return levels.merge_by_energy_order(primary.identified, other.identified)


def _relabel(other, merged):
    """The second folder, keyed so a level of the first can find its twin."""
    other.twin = {}
    for one in merged:
        if one.matched:
            other.twin[(one.primary.kset, one.primary.j,
                        one.primary.rcg)] = one.secondary.rcg
    return other


def _report_energy_order(merged):
    worst = levels.energy_order_report(merged)[:5]
    for (kset, j), one in worst:
        print('   set %d J=%-4s worst dE = %+.4f cm-1'
              % (kset, ('%.1f' % j), one.delta))
    missing = levels.unmatched(merged)
    if missing:
        print('   %d levels of the first run found no twin' % len(missing))


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('printout: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
