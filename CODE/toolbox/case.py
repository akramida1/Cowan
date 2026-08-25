"""One calculation folder, with everything in it tied together.

A Cowan case is six or seven files that only make sense against each
other: IN36 says what the configurations are, ING11 says which subshell
slot is which, OUTG11 carries the eigenvectors, genealogy.txt the
coupling chains, and the RCE files the fitted energies.  Every script in
this package needs the same joins between them, so they are made once,
here, and `printout.py` and `conv_out.py` ask this object for what they
need.

Nothing is read that is not asked for: the OUTG11 pass is filtered to
the couplings the caller named, which is what keeps a Pr III run to a
few seconds instead of a walk through 122 MB per coupling.
"""

from dataclasses import dataclass, field as _field
from pathlib import Path

from . import designate, files, genealogy, in36, ing11, levels, outg11
from . import parnames, rce


class CaseError(Exception):
    """A folder that does not hold a complete, consistent run."""


@dataclass
class Case:
    """One folder, read once.

    `parent` is the folder to take IN36 from, for a sub-folder that
    carries only the files RCG and RCE rewrite -- the `JK/` arrangement
    the author uses to get a second coupling out of the same fit.  It
    need not be given: a folder with no IN36 of its own falls back to
    the directory above it (`_in36_folder`), which is that same
    arrangement seen from inside the sub-folder.
    """

    folder: object
    couplings: tuple = ('LS',)
    parent: object = None
    drop: tuple = ()                    # (n, l) to leave out where closed

    deck: object = None
    configurations: object = None       # in36.In36
    core: tuple = ()                    # the closed shells IN36 shares
    shells: dict = _field(default_factory=dict)     # (kset, index) -> slots
    blocks: dict = _field(default_factory=dict)     # (kset, J) -> JBlock
    genealogy: object = None
    csets: list = _field(default_factory=list)      # rce.Cset
    starting: dict = _field(default_factory=dict)   # OUTGINE's parameters,
                                                    # where RCE has not run
    identified: list = _field(default_factory=list)
    unlabelled: set = _field(default_factory=set)   # (coupling, kset, J) that
                                                    # RCG printed in LS only
    unmatched: list = _field(default_factory=list)  # RCE levels above the
                                                    # deck's E0MAX that RCG
                                                    # never printed
    _named: dict = _field(default_factory=dict, repr=False)
    _designations: dict = _field(default_factory=dict, repr=False)

    # ----------------------------------------------------------------
    # loading
    # ----------------------------------------------------------------

    @classmethod
    def read(cls, folder, couplings=('LS',), parent=None, with_rce=True,
             drop=(), use_emax=False):
        case = cls(folder=folder, couplings=tuple(couplings), parent=parent,
                   drop=tuple(drop))
        case.deck = ing11.read_ing11(folder=folder)
        case.configurations = in36.read_in36(folder=_in36_folder(folder,
                                                                parent))
        case.core = in36.common_closed_shells(case.configurations.configurations)
        by_label = case.configurations.by_label()
        # The order the slots hand out their `n` values: Cowan's own
        # left to right for a deck RCN2 wrote, the permutation the cards
        # agree on for a reordered one (`ing11.slot_order`).
        order = ing11.slot_order(case.deck.configurations, by_label, case.core)
        for card in case.deck:
            try:
                configuration = by_label[card.label]
            except KeyError:
                raise CaseError(
                    'the configuration %r is on an ING11 card of %s but not '
                    'in its IN36' % (card.label, folder))
            case.shells[(card.kset, card.index)] = ing11.match_shells(
                card, configuration, case.core, order)

        case._check_couplings()
        # LS is always read whatever was asked for: it is the scheme the
        # RCE-to-RCG identification runs on (`levels.identify`), because
        # it is the one both programs write eigenvector components in.
        wanted = list(case.couplings)
        if designate.LS not in wanted:
            wanted.append(designate.LS)
        case.blocks = levels.by_block(outg11.read_outg11(folder=folder),
                                      wanted)
        case.genealogy = genealogy.read_genealogy(folder=folder)
        if not with_rce:
            # The deck's own parameters, under the names RCG punched for
            # RCE.  Without them a folder RCE has never run in can be
            # reported on for its levels but not for its parameters,
            # which are in ING11 all along but nameless there.
            case.starting = rce.read_outgine(folder=folder)
        if with_rce:
            case.csets = rce.read_rceout(folder=folder)
            case._refine()
            emax = case.deck.emax if use_emax else (0.0, 0.0)
            for cset in case.csets:
                limit = emax[cset.index - 1] if cset.index <= 2 else 0.0
                case.identified.extend(
                    levels.identify(cset, case.blocks, emax=limit or None,
                                    skipped=case.unmatched))
        return case

    def _refine(self):
        """Take the energies and the Lande g factors from LEVELS1.

        RCEOUT prints the fitted energy as F10.3 and no g at all; LEVELS1
        prints the same levels in the same order with F11.4 and the g
        beside them (`FOR/RCE20KD.F:2033`).  RCEOUT is still the file the
        identification runs on, because it is the one that lists the
        eigenvector components.

        LEVELS1 repeats each parity set once per coupling scheme, so the
        sections are matched to the parity sets by length and then
        checked level by level on J; a disagreement leaves the RCEOUT
        values alone rather than writing the wrong energies.
        """
        try:
            sections = rce.read_levels(folder=self.folder)
        except (OSError, IOError, rce.RceError):
            return
        at = 0
        for cset in self.csets:
            while (at < len(sections)
                   and len(sections[at].levels) != len(cset.levels)):
                at += 1
            if at >= len(sections):
                return
            pairs = list(zip(cset.levels, sections[at].levels))
            if any(abs(a.j - b.j) > 1e-6 for a, b in pairs):
                return
            names = sections[at].configurations
            for a, b in pairs:
                a.e_obs, a.e_calc, a.g = b.e_obs, b.e_calc, b.g
                a.conf0 = self._full_label(b.conf0, names)
            at += 1

    def _full_label(self, truncated, names):
        """The full ING11 label behind a LEVELS dominant-configuration field.

        RCE writes that field as A6 (`FORMAT` at `FOR/RCE20KD.F:1988`)
        while the label it names is the A8 of ING11, so a name longer
        than six characters comes through cut.  The section's own header
        lists the full names of exactly the configurations that can
        appear in it (`FOR/RCE20KD.F:1305`), so the cut is undone by
        matching against those and nothing wider.

        Two names of one parity set agreeing in their first six
        characters would make the field genuinely ambiguous.  That stops
        the run by name rather than picking one: a wrong configuration
        here silently regroups a purity table.
        """
        if not truncated:
            return ''
        found = [n for n in names if n[:6].strip() == truncated]
        if len(found) == 1:
            return found[0]
        if not found:
            # Not in the header at all: keep what RCE wrote.  A deck can
            # be read without its LEVELS header being complete, and a
            # six-character label is already whole.
            return truncated
        raise CaseError(
            'the LEVELS dominant-configuration field %r of %s could be any '
            'of %s: RCE prints it as A6 and these labels do not differ '
            'within six characters, so the level cannot be assigned to a '
            'configuration' % (truncated, self.folder, ', '.join(found)))

    def _check_couplings(self):
        """Every coupling asked for must be one the deck asked RCG for.

        The plan is explicit that a missing coupling is a hard error and
        not something to work around: the user names what they want, and
        a silently short report is worse than none.
        """
        asked = self.deck.main.couplings
        missing = [c for c in self.couplings if c not in asked]
        if missing:
            raise CaseError(
                '%s does not have %s.  Its ING11 main card asks RCG for %s; '
                'the couplings are columns 31-39 of that card, and a digit '
                'of 0 turns one ON.  True jj is column 38 and is the one '
                'exception -- it needs KCPLD(8)=8.'
                % (self.folder, ', '.join(missing),
                   ', '.join(asked) or 'nothing'))

    # ----------------------------------------------------------------
    # what the reports ask for
    # ----------------------------------------------------------------

    def card(self, kset, index):
        return self.deck.set(kset)[index - 1]

    def designation(self, coupling, kset, j, state):
        """The ASD columns of one basis state of one J block."""
        key = (coupling, kset, j, state)
        try:
            return self._designations[key]
        except KeyError:
            pass
        outer = None
        if coupling in designate.K_SCHEMES:
            # Only these need the eigenvector row: K is theirs alone, and
            # no genealogy table carries it.
            try:
                section = self.blocks[(kset, j)].section(coupling)
            except outg11.Outg11Error:
                if self.blocks[(kset, j)].nlevels != 1:
                    raise
                # An OUTG11 written before RCG printed the sections
                # of a one-state block has this one in LS and in
                # nothing else, so K was never written.  The level is
                # that one basis state whatever the coupling, so it is
                # named in LS rather than guessed at; `unlabelled`
                # counts them so the report can say so.
                self.unlabelled.add((coupling, kset, j))
                coupling = designate.LS
            else:
                outer = designate.read_outer(coupling,
                                             section.rows[state - 1].desc)
        gb = self.genealogy.block(kset, j)
        chain = gb.state(designate.TABLE_OF[coupling], state)
        found = designate.designate(coupling, chain,
                                    self.shells[(kset, chain.cfg)],
                                    outer=outer, drop=self.drop)
        self._designations[key] = found
        return found

    def amplitudes(self, coupling, kset, j, eigenvector):
        """(|amplitude|, basis state) of one eigenvector, largest first.

        RCG used to print no section for a block of one basis state
        in any coupling but LS: there is nothing to transform, so TMX
        was never formed and the print loop never reached.  It prints
        them all now -- true jj first, with the LSJJRD fix, and then
        couplings 3-7 -- but an OUTG11 written before that does not
        have them, which is why this fallback stays.  The row is
        recovered rather than dropped -- a
        one-dimensional block has amplitude 1 by construction, and
        `genealogy.txt` prints its chain whatever RCG left out.  A block
        of more than one state with no section is a real gap and is
        reported as one.
        """
        block = self.blocks[(kset, j)]
        try:
            section = block.section(coupling)
        except outg11.Outg11Error:
            if block.nlevels != 1:
                raise
            return [(1.0, 1)]
        weights = [(abs(row.amplitudes[eigenvector - 1]), row.number)
                   for row in section.rows
                   if eigenvector <= len(row.amplitudes)]
        weights.sort(reverse=True)
        return weights

    def parameters(self, kset):
        """The parameters of one parity set, under readable names."""
        if kset not in self._named:
            self._named[kset] = parnames.resolve(
                self.parameter_source(kset), self.deck.set(kset),
                {c.label: self.shells[(c.kset, c.index)] for c in self.deck})
        return self._named[kset]

    def parameter_source(self, kset):
        """The parameter table the reports run on, in RCE's own order.

        RCEOUT where RCE has run, and OUTGINE where it has not: the two
        list the same parameters in the same order, and differ in that
        OUTGINE's values are the deck's rather than a fit's and carry no
        standard deviations.
        """
        if self.csets:
            return self.cset(kset).parameters
        try:
            return self.starting[kset]
        except KeyError:
            raise CaseError('%s has no parameters for parity set %d: RCEOUT '
                            'was not read and OUTGINE does not have it'
                            % (self.folder, kset))

    def cset(self, kset):
        for cset in self.csets:
            if cset.index == kset:
                return cset
        raise CaseError('%s has no RCE parity set %d' % (self.folder, kset))

    def scale_factors(self, kset, index):
        """The HR factors a configuration's parameters were scaled by.

        `IFACT`, columns 73-80 of the parameter card, as percentages, in
        the order Fk(i,i), zeta, Fk(i,j), Gk(i,j) (`rcg11kd.f:10837`);
        a blank or 99 means 100.  They are needed to divide a
        Hartree-Fock value back out of the scaled one RCG punched.
        """
        return _factors(self.deck.parameters[kset - 1][index - 1])

    def ci_factor(self, kset, labels):
        """The one HR factor a configuration-interaction card carries.

        A CI card has a single scale factor, `IFACT(4)` -> `FACT5`
        (`rcg11kd.f:10905-10917`), not the four of a parameter card.
        """
        for card in self.deck.ci[kset - 1]:
            if set(card.labels) == set(labels):
                return _factors(card)[3]
        return 1.0


def _in36_folder(folder, parent):
    """Where to read IN36 from: the run's own folder, or the one above.

    RCG and RCE need no IN36 -- it is RCN's input -- so a sub-folder
    holding a second coupling of the same fit does not have one, and the
    author runs `printout` in it as readily as in the folder above.  An
    explicit `parent` still wins, because a caller that names one has
    decided; otherwise the directory above is tried, and only if it has
    no IN36 either does the missing-file error come back naming the
    run's own folder, which is where a reader would look first.
    """
    if parent is not None:
        return parent
    if files.resolve(files.IN36, folder, required=False) is not None:
        return folder
    above = Path(folder).resolve().parent
    if files.resolve(files.IN36, above, required=False) is not None:
        return above
    return folder


def _factors(card):
    """A card's four HR percentages as fractions; blank or 99 means 100."""
    return tuple(1.0 if f in (0, 99) else 0.01 * f for f in card.factors)
