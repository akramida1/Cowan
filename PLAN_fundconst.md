# Plan — factor all fundamental constants into `FOR/fundconst.f`

Status: in progress. Written 2026-08-20.  Sections 0c, 1, 2, 3, 4, 5 and 6
are done, and the three rebuilt executables are committed (see "Promotion"
after Section 6); Section 7 (the RCG conversion-coefficient routines) is next.  The
constants pass over the Fortran is finished except for RCE, which Section 8
expects to have nothing in it.

Scope: replace every hard-coded physical constant in the Fortran with named
parameters from a single include file, update those values to CODATA 2022,
and add exact E1/M1/E2 conversion coefficients to RCG.

---

## 0. The reported gf/gA discrepancy — resolved, and it is not a constant

`WORK/gA_gf_check.txt` shows ±1e-3 scatter of random sign between the printed
`GA` and the `GA` recomputed from the printed `LOG GF`. That scatter is
**quantization of the 3-decimal `LOG GF` print field**, not an arithmetic error.

A half-unit in the last printed decimal of a base-10 logarithm is
`ln(10) × 0.0005 = 1.151e-3` in relative gf. Recomputing the exact log gf
implied by the printed `GA` and `dE`, using the modern
`8·pi^2·alpha^2·a0·c = 0.66702518`:

| sigma (cm-1) | log gf printed | log gf implied | delta | half-ulp = 0.0005 |
|---|---|---|---|---|
| 418430 | -1.485 | -1.48521 | +0.00021 | inside |
| 411695 | -2.049 | -2.04940 | +0.00040 | inside |
| 408365 | -2.745 | -2.74521 | +0.00021 | inside |
| 397156 | -2.285 | -2.28511 | +0.00011 | inside |
| 396287 | -3.276 | -3.27563 | -0.00037 | inside |
| 384932 | -4.247 | -4.24688 | -0.00012 | inside |
| 368063 | -4.965 | -4.96458 | -0.00042 | inside |

All seven lines are inside ±0.0005. The printed `LOG GF` is the correctly
rounded value of the gf that is exactly consistent with the printed `GA`.
Every observation matches: random sign, ~0.07% rms, ~0 mean, second decimal
affected.

Consequences for this plan:

- Updating constants will **not** remove this scatter. There is nothing to fix
  in the arithmetic.
- **Decision (author, 2026-08-20): no numerical correction is to be made.**
  The `LOG GF` print field is *not* to be widened — that would move `OUTG11`
  columns and break the Python and Perl parsers, to fix something that is not
  broken. What is required instead is a **prominent, permanent reminder** in
  three places that in an RCG transition list **`gA` is the primary reference
  quantity, not `log gf`**. See Section 0c.
- The constants work below is still worth doing on its own merits.

## 0c. The `gA`-is-primary reminder — **DONE 2026-08-20**

This has been investigated before and the result keeps being forgotten, so it
must leave a permanent footprint in the source and the documentation rather
than in anyone's memory.

**The point to record**, including the reason the wrong intuition is so
persistent:

> In the transition lists printed by RCG, **`gA` is the primary reference
> quantity; `log gf` is not.** `log gf` is printed to three decimals, and a
> half-unit in the last decimal of a *base-10 logarithm* is
> `ln(10) x 0.0005 = 1.15e-3`, i.e. **0.115 % in gf** — already in the second
> printed decimal when converted back. Recovering `S`, `f` or `A` from
> `log gf` therefore injects ~0.1 % of random-sign scatter that is not present
> in `gA`.
>
> The naive expectation is the opposite — that `gf` should be the more
> accurate of the two, because it is obtained from `S` with a *lower* power of
> the transition energy (sigma^1 rather than sigma^3), so energy rounding
> propagates into it less strongly. That reasoning is correct as far as it
> goes, and it is beside the point: the dominant error is not in the energy,
> it is contributed by the `log10` at print time, and that step is routinely
> overlooked. `gA` is printed in floating form and carries its own digits;
> `log gf` does not.

Placement — all three are required, and the wording is to be the same
argument in each so that whichever one is found first is sufficient:

1. **`README.md`** — in the section describing RCG output, prominent.
2. **Top of `FOR/rcg11kd.f`** — in the header comment block, positioned so it
   is visible without scrolling into the change list.
3. **Top of `SUBROUTINE SPECTR` (`FOR/rcg11kd.f:12889`)** — the routine that
   actually writes the transition list, immediately after the `SUBROUTINE`
   statement.

Comment-only; bit-exact by construction. **Done:** 40 comment lines
added to `FOR/rcg11kd.f` (header block before the change list, and after the
`SUBROUTINE SPECTR` statement, now at line 12910) and 20 lines to `README.md`
under item 3 of ADDITIONS TO THE MANUAL. `git diff --stat` = 60 insertions,
0 deletions; every added Fortran line begins in column 1 with `C` and none
exceeds column 72; `gfortran -fsyntax-only` with the full flag set passes.

## 0a. Two genuine defects found while surveying

1. **`FOR/rcg11kd.f:13312`** — `GA=0.6698D0*PSQ*ERAS*(TEMP**3)`. A dropped
   digit; all seven other sites write `0.66698D0`. It is +0.42% and drives the
   `TEXCIT`-weighted `D0` line-selection threshold in `SPECTR`, so it silently
   changes which lines survive the cutoff. Not a printed quantity.
2. **`FOR/RCN36K.F` uses two different fine-structure constants.** `zeta1`
   uses `2.66246e-5` = alpha²/2 from a pre-1969 alpha (−4.0e-5 relative),
   while `sm`/`sn`/`vk` use `2.921830` = (alpha²/2)·Rinf evaluated with the
   *modern* alpha, exact to all eight digits. Half of it was updated at some
   point and half was not.

## 0b. The existing structure is right; only the numbers are stale

Verified algebraically against the Grant / Martin–Wiese formulas in
`WORK/A-S-f-conversion_formulas.docx`:

- `PSQ = 3.0375D-6 * PMAX**2` is `2/(3·2Rinf)`, the E1 gf coefficient. Correct
  form; the value is 1.9e-5 low.
- `IF (II.EQ.0) PSQ=PSQ*(0.5D0/137.04D0)**2` is `(alpha/2)²`. That is exactly
  the M1/E1 ratio: f(M1)/f(E1) = alpha²/4. Correct form.
- `IF (II.EQ.2) PSQ=PSQ*((3.14159D0*0.529D0)**2)/5.D0` combined with
  `GF=GF/(ERAS2**2)` is `pi²·a0²·sigma²/5`. That is exactly the E2/E1 ratio.
  Correct form.
- `GA = 0.66698D0*GF*sigma**2` is `8·pi²·alpha²·a0·c`, whose correct value is
  0.66702518. Correct form, 6.8e-5 low.

So this is a pure constants job. No formula is wrong.

---

## 1. Inventory

57 literal sites in 19 routines out of 178 across the package.

| File | Sites | Routines |
|---|---|---|
| `FOR/RCN36K.F` (RCN) | 22 | `rcn`, `scheq`, `phshift`, `sli1`, `zeta1`, `outpt`, `zetabw`, `sm`, `sn`, `vk` |
| `FOR/RCN2K.F` (RCN2) | 14 | `G5INP`, `SLI2`, `ZETA2` |
| `FOR/rcg11kd.f` (RCG) | 21 | `RCOEFF`, `SPECTR`, `ENERGY`, `BORN`, `LVDIST`, `WNDIST` |
| `FOR/RCE20KD.F` (RCE) | **0** | pure least-squares fitter; the only decimals are in comments |
| `LSJJ*.F`, `RCEJJ2E.F` | 0 | — |

Distinct constants and their error against CODATA 2022:

| Literal | Meaning | Rel. error |
|---|---|---|
| `109737.31` `109.73731` `109737.3` `109737.0` | Rinf | 5e-8 … 3e-6 |
| `8065.479` `8065479.` | eV → cm-1 | −8.0e-6 |
| `137.036` `137.0360` `137.04` | 1/alpha | −7e-9, −2.6e-5 |
| `274.072` | 2/alpha | +3e-9 |
| `2.66246e-5` `1.33123e-5` | alpha²/2, alpha²/4 | −4.0e-5 |
| `2.921830` `5.843660` | (alpha²/2)·Rinf, ×2 | 0 — already exact |
| `3.0375D-6` | 2/(3·2Rinf), E1 gf | −1.9e-5 |
| `0.66698` (×7) | 8·pi²·alpha²·a0·c | −6.8e-5 |
| `0.6698` (×1) | same | **+4.2e-3** (defect 0a.1) |
| `0.529` | a0 in 1e-8 cm | −3.3e-4 |
| `0.0136058` `13.605` | Ry in eV | −5e-6 … −5.1e-5 |
| `3.14159` `3.141592654` `3.141592653590` | pi | 8.5e-7 … 0 |

---

## 2. `FOR/fundconst.f` — design

Fixed form, `PARAMETER` statements only — no `COMMON`, no executable
statements — so it can be included anywhere after `IMPLICIT`. `INCLUDE` is
already the house idiom (80 sites in `rcg11kd.f`, 37 in `RCE20KD.F`).

All four programs are `IMPLICIT REAL*8 (A-H,O-Z)` throughout, so one
double-precision file serves every one of them.

**Contents** — the CODATA 2022 base set in SI, plus the CGS/atomic duplicates
the code actually uses, plus derived conversion coefficients. Each line
carries its CODATA value and uncertainty in a trailing comment, and values
that are exact by definition (`h`, `e`, `c`, `k`, `NA`) are marked as such.

- Base: `CLIGHT` (cm/s and m/s), `PLANCK`, `ECHARG`, `EMASS`, `EPS0`, `AMU0`,
  `BOHRMG`, `ABOHR` (cm and m), `ALPHFS`, `RINF` (cm-1), `HARTJ`.
- Derived: `ALPHF2` = alpha², `AUCM` = 2·Rinf, `EV2CM`, `RYEV`, `PICON`,
  `GFE1`, `GACON`.

**Name collisions — checked against every identifier in all four programs:**

- `ALPHA` is already an ordinary variable, 19 uses. Use `ALPHFS`.
- `PI` is already an ordinary variable, 9 uses. Use `PICON`.
- `CLIGHT`, `RINF`, `AUCM`, `EV2CM`, `ABOHR`, `RYEV`, `GFCON`, `GACON`: zero
  hits anywhere. Safe.

Re-run the collision check before each file's substitution. A `PARAMETER` name
that shadows a local is a silent wrong-answer bug, not a compile error, when
the local is only read.

---

## 3. The acceptance-standard problem, and the two-stage protocol

CLAUDE.md §1 makes byte-identical output the acceptance test. This change
deliberately changes numbers, so §1 cannot apply directly — and if the
refactor and the value update land in one pass, a plumbing mistake (wrong
parameter substituted, precision change, name shadowing) becomes
indistinguishable from an intended constant change. There would be no oracle
at all.

Every file below is therefore done in two stages:

- **Stage A** — introduce the parameter *carrying the legacy value, bit for
  bit as written today*, and substitute mechanically. Output must be
  byte-identical. §1 applies in full. This proves the plumbing.
- **Stage B** — change only the value, in `fundconst.f`. Every diff that then
  appears is by construction attributable to that constant.

Stage A is where the bugs get caught. Do not skip it to save a build.

**Stage A hazard specific to `RCN36K.F`:** its literals are written without a
`d0` suffix (`137.036`, `2.66246e-5`, `2.921830`). Those are *default real*
literals — evaluated in single precision and promoted at use. Substituting a
double-precision `PARAMETER` changes the result **even when the decimal digits
match**. Stage A will not be bit-exact for those sites unless the parameter is
built to reproduce the single-precision rounding. Enumerate these ~12 sites up
front and decide each one explicitly, rather than letting them hide inside the
Stage B diff.

---

## 4. Section-by-section work

Order is top-to-bottom through the pipeline, as requested: RCN → RCN2 → RCG →
RCE. Each section is independently reviewable and independently revertable.

### Section 1 — Numeric comparison tool (prerequisite) — **DONE 2026-08-20**

`diff` stops being an oracle the moment Stage B runs. Build a tolerance-based
comparator first.

- `CODE/cmp_num.py` + `CODE/cmp_num.bat` (lowercase, per CLAUDE.md §3), over
  `CODE/toolbox/cmpnum.py`, following the script-over-toolbox split every
  other utility here uses.
- Tokenizes both files into a text skeleton plus a numeric vector. Requires
  the text skeletons to match exactly — that catches structural change, which
  is never acceptable. Then compares the numeric vectors with relative and
  absolute tolerances (a pair passes on *either*, so that fields printed to a
  fixed number of decimals near zero are not judged on a relative deviation
  they do not have).
- Reports: count compared, max relative deviation with file/line/column, rms,
  a decade histogram, and the N worst offenders. Every pair enters the
  histogram whether or not it passed — the shape of the distribution is the
  deliverable.
- Exit codes are three-valued and the distinction is load-bearing:
  **0** within tolerance, **1** numbers moved and the text is intact,
  **2** structural mismatch or a missing file. 1 says the change did what it
  was meant to; 2 says something moved that was never allowed to.
- Accepts two files or two folders. `--rtol --atol --top --ignore
  --no-default-ignores --d-exponents --quiet --only`.
- The three CLAUDE.md §1 timing lines (`TIME=`, `FINISHED … AT`, `DIAGTIME`)
  plus RCN's `t(run)=` / `t(job)=` are ignored by default.
- 32 tests in `CODE/tests/test_cmpnum.py`, all fast; they run against the real
  Sn VIII outputs in `WORK/`, so no `slow` marker was needed.

Four tokenizer decisions that are worth knowing before reading a report:

1. **A label is matched first and stepped over.** `Sn7+`, `d65p` and the `d7`
   of `4d7` stay in the skeleton. A lookbehind cannot do this: reject a digit
   after a letter and the *next* digit of `d65p` still comes back as `5`;
   also reject a digit after a digit and the sign of the second of two packed
   fields is silently dropped, so `418.430-238.989` compares as `+238.989`.
   Both were caught by the tests, and both were real.
2. **`E` is the only exponent letter; `D` is not**, unless `--d-exponents`.
   No output file the package writes has a `D` exponent — checked across
   OUTG11, OUT36, OUT2, OUTGINE, TAPE2E — while `4d7` is everywhere, and
   reading it as an exponent makes it 4e7.
3. **A bare-sign exponent (`1.234-105`) is read as two numbers.** Ambiguous by
   construction; resolved towards the case that actually occurs here.
4. **Blank runs collapse in the skeleton**, because a number that gains a
   digit eats a blank of its fixed-width padding, and that must not read as a
   structural change on top of the numeric one.

Fully packed fields with no separator at all would be read as one long number.
No FORMAT in use produces one — a test asserts that across all five outputs —
and `long_tokens` in the report counts them if one ever appears.

**`TAPE2N` is Fortran unformatted and is not tokenized.** It gets an md5
answer — identical, or differing by an amount the tool will not pretend to
quantify — and is never silently skipped. See the note in Section 4.

This tool is also what makes the RCN result reviewable at all, since RCN's
output becomes the new benchmark for everything downstream.

### Section 2 — Baselines — **DONE 2026-08-20**

Run, from the current unmodified build, and keep md5s plus the raw files:

- `WORK/` (Sn VIII) — RCN + RCN2 + RCG; small and fast, and the only case that
  exercises `sm`/`sn`/`vk`/`phshift`.
- `WORK/Pr3/` — the §1 regression case.
- `WORK/Forbid` — the only M1/E2 deck. Required for Section 6.

Per the `tape2n-differs-across-platforms` note, `TAPE2N` carries platform
rounding noise. Baseline on this machine and do not chase it.

**Done.** Everything is in **`F:\COWAN\fundconst_base\`**, outside the
repository — outputs are not committed and `Pr3/OUTG11` alone is 85 MB.
`README.txt` there is the full account; the short form:

- `bin/` — the four programs rebuilt from the sources at commit `2ed5631`
  with the CLAUDE.md §2 flag sets (gfortran 13.2.0 MinGW-W64, Strawberry),
  the only warning being the expected `-fno-automatic overwrites
  -frecursive`. Beside them: `SENIOR`, `LSJJ.DAT`, `ING11.CFP`, the four
  DLLs, and `FOR072`/`FOR073`/`FOR074`. The CFP decks are not in the
  repository and are machine-built, so they were made **by this build**,
  following `make_cfp.bat` by hand.
- `SnVIII/`, `Pr3/`, `Forbid/` — one directory per case, `IN2` + `IN36`
  copied in and run by `run_case.sh` (RCN, RCN2, RCG; `Forbid` also RCE).
  Console output kept per program. `MD5SUMS.txt` covers inputs and outputs
  of all three.

Three things were established rather than assumed:

1. **The baseline build is a valid oracle.** `SnVIII/` against the committed
   `WORK/`: `OUT2`, `ING11`, `OUTGINE`, `TAPE2E`, `TAPE2N` md5-identical;
   `OUT36` and `OUTG11` differ **only** in the §1 timing lines. `cmp_num.py`
   compares 886 numbers in `OUT36` and 53148 in `OUTG11` at max relative
   deviation exactly 0, and answers for `TAPE2N` with `md5 identical`. That
   is also the first exercise of the Section 1 tool on real files of both
   kinds. The executables themselves are *not* byte-identical to the
   committed `CODE/*.exe` — a gfortran build is not reproducible to the byte
   — but their output is, which is what the standard is about.
2. **`Forbid/` really does carry the forbidden lines**: 89 M1 and 131 E2,
   exactly as `WORK/Forbid/README.txt` documents. Section 6 has something to
   watch.
3. **`Pr3/` is the publication deck, not the smaller regression deck** of §1
   — its own README says so, and adds that this sequence's `TAPE2E` is
   byte-identical to the regression run's. RCG wrote the 85 MB `OUTG11` that
   README predicts. RCE was not run: Section 2 needs RCN/RCN2/RCG only.

Re-running all three cases costs under fifteen seconds (Sn VIII ~1 s,
Pr3 ~11 s, Forbid ~1 s), so no Stage A check has an excuse to skip one.

### Section 3 — `fundconst.f` written, nothing yet included — **DONE 2026-08-20**

File authored with legacy values only, one commented CODATA 2022 value per
line. Compiles as part of a throwaway include test. No behaviour change.

Reviewable on its own, and the whole point of the exercise: this one file
becomes the single place a physicist has to look.

**Done.** `FOR/fundconst.f`, 323 lines, untracked-new; nothing else touched.
Section A is 16 CODATA 2022 constants that nothing uses yet (the SI base set
plus the CGS/atomic duplicates, `ALPHF2`, `RINFCD`, `AUCM`); Section B is the
23 names that will replace literals, each carrying its **legacy** value, with
the CODATA value, the relative error, the list of call sites, and — where the
name exists only to reproduce a coarser rounding — the expression it collapses
to at Stage B.

- **Values came from NIST, not from memory** (open question 2 is closed for
  every value in the file): `physics.nist.gov/cuu/Constants/Table/allascii.txt`,
  2022 adjustment. Derived quantities were computed from those digits. Three
  §0b/§7 numbers are confirmed to the digit: `2/(3·AUCM)` = 3.0375568352755e-6,
  `8·pi²·alpha²·a0·c` = 0.667025176604160, `(alpha²/2)·Rinf` = 2.9218303467490
  — the last confirming §1's "already exact" for `2.921830`.
- **Every name is declared `DOUBLE PRECISION` explicitly**, not left to the
  host routine's `IMPLICIT REAL*8 (A-H,O-Z)`. Not every unit of `RCN36K.F`
  carries that card, and `KBOLTZ` starts with K, so under `IMPLICIT INTEGER*4
  (I-N)` it would otherwise be an integer. Verified both ways.
- **Collision check re-run** over `FOR/*.F`, `FOR/*.f` and `FOR/LSJJ/`, word
  boundaries, case-insensitive: **zero hits for all 39 names.**
- **Include test.** A throwaway program prints all 23 legacy values beside
  their CODATA targets and checks the internal identities: `AUCM-2*RINFCD`
  and `CLIGHT-100*CLIGHM` and `ABOHR-100*ABOHRM` are exactly zero, and
  `ALPHF2*(1/alpha)²-1` is -8.6e-12. It also includes the file into a routine
  with **no** `IMPLICIT` card and into a lowercase `implicit real*8` routine.
  Beyond that, `fundconst.f` was inserted after the first `IMPLICIT` of each
  of **all four real sources** and `-fsyntax-only` run with the full §2 flag
  set: four clean passes, the `-fno-automatic overwrites -frecursive` warning
  being the only output. Those temporary copies were deleted.
- Fixed form verified mechanically: every line is a `C` comment in column 1 or
  starts in column 7, and none exceeds column 72.

**Two corrections to the §1 inventory**, found while writing the file. Neither
changes the design; both change what Sections 4 and 6 have to do.

1. **The single-precision hazard is not ~12 sites, it is all 22 of
   `RCN36K.F` — plus one in `rcg11kd.f`.** `RCN36K.F` uses a `d0` suffix on
   *none* of its constant literals, so every one of them is a default real,
   rounded to single precision and promoted only where it meets a double.
   Substituting a double `PARAMETER` moves each by ~3e-8 relative — larger
   than the `RINF`, `ALPHIN`, `TWOOAL`, `ZETCON`, `ZETCN2`, `PICON` and
   `FOURO3` updates it is supposed to be separable from. The overlooked
   `rcg11kd.f` site is line 11046, `ERAS=ERAS*109737.0/UENRGY`, the one
   constant in that file written without `D0`. All 23 are enumerated by line
   number in the file's header warning. `RCN2K.F` is unaffected.
2. **Live-site counts differ from the table**, because two of the counted
   sites are inside commented-out code: `rcg11kd.f:11968` (`3.14159D0`, in the
   commented-out skewed-Gaussian plot in `LVDIST`) and `RCN2K.F:2318`
   (`PI=3.141592654D0`). So **`LVDIST` has no live site at all** — RCG's 6
   routines are really 5: `ENERGY`, `SPECTR`, `WNDIST`, `BORN`, `RCOEFF` —
   and RCN2's pi is not a site either. Recount of live literals: RCN36K 22,
   RCN2K 13, RCG 26, RCE 0 = **61**, against the table's 57.

No parameter carrying a single-precision-rounded value was defined, on
purpose: §4 requires each of those 23 sites to be decided explicitly rather
than absorbed into the constants diff.

### Section 4 — RCN (`RCN36K.F`) — **DONE 2026-08-20**

10 routines, 22 sites, including every single-precision-literal hazard and
the two-alphas defect (0a.2). The riskiest file, and the one whose output
feeds everything else.

**Three stages, not two.** §3's Stage A/Stage B split cannot absorb the
default-real problem: promoting a literal to a double `PARAMETER` moves the
result by ~5e-8 even with identical digits, which is *larger* than the
`RINF`, `ALPHIN`, `TWOOAL`, `ZETCON` and `ZETCN2` updates it was meant to be
separable from. So the precision change was given a stage of its own, and
every one of the 22 sites was decided the same way and for the same reason:
there is no site among them where the single rounding is doing anything but
losing digits on the way into double arithmetic. It is an artefact of the
source predating `d0` suffixes, not a choice. The reasoning is recorded in
the `RCN36K.F` header (Vers. 36.4.2) rather than only here.

| Stage | What changed | Result |
|---|---|---|
| A | literal → `REAL` parameter with the same decimal string | **byte-identical** |
| B1 | those → the Section B `DOUBLE PRECISION` names | precision only |
| B2 | Section B values → CODATA 2022 | values only |
| B3 | `RINFL`→`RINF`, `PICONB`→`PICON`, now equal by definition | **byte-identical** |

Before any of it, the unmodified sources were rebuilt into a scratch tree and
shown to reproduce the Section 2 baseline exactly — `TAPE2N` md5-identical and
`OUT36` differing only in `t(run)`/`t(job)` — on all three cases. Without
that the Stage A result would have proved nothing about the edit.

Stage A and Stage B3 are byte-identical on `SnVIII`, `Pr3` and `Forbid`
(`TAPE2N` md5-identical, `OUT36` identical apart from the timing lines), so
CLAUDE.md §1 applies to them in full.

**Stage B1 — the precision promotion alone.**

| case | numbers | max rel | rms | shape |
|---|---|---|---|---|
| SnVIII | 886 | 5.14e-08 | 5.2e-09 | 859 exact, the rest at 1e-8..1e-7 |
| Forbid | 423 | 4.79e-08 | 5.1e-09 | 408 exact, the rest at 1e-8..1e-7 |
| Pr3 | 22639 | 6.39e-02 | 7.1e-04 | 21176 exact, bulk at 1e-9..1e-7 |

Single-precision epsilon and nothing else. Pr III's three outliers are
`deltn(118..120)`, the SCF convergence residuals — 2e-6 quantities whose
own scale *is* the noise floor, so a 6% move in them is the expected
behaviour of an iterative solver under a 5e-8 perturbation, not a finding.

**Stage B2 — the values alone.**

| case | numbers | max rel | rms | bulk |
|---|---|---|---|---|
| SnVIII | 886 | 1.15e-04 | 1.3e-05 | 1e-5..1e-4 |
| Forbid | 423 | 1.14e-04 | 1.3e-05 | 1e-5..1e-4 |
| Pr3 | 22639 | 6.39e-02 | 7.1e-04 | 1e-5..1e-4 |

The two effects are cleanly separated: B2 is ~2000x B1 in the two small
cases, so nothing in the constants diff is precision noise in disguise.

**What actually moved, and it is not the Rydberg constant.** `RINF` is
5.2e-8 stale and `PICON` 6.6e-14; the change is dominated by defect 0a.2,
the two different fine-structure constants. `zeta1` built the spin-orbit
integral from `2.66246e-5` — alpha²/2 from a pre-1969 alpha, 4.05e-5 low —
while `zetabw`, `sm`, `sn` and `vk` used `2.921830` = (alpha²/2)·Rinf from
the *modern* alpha, correct to every digit it carried. Both now come from
the same `ALPHFS`. The 1e-5..1e-4 bulk above is that 4.05e-5 and its
propagation through the SCF; the 1e-4 tail is derived quantities carrying
it with a small amplification.

**`TAPE2N`.** `CODE/toolbox/tape2n.py` was written, mirroring the one
`write (2)` in `subroutine outpt`, and `cmp_num.py` now reads TAPE2N by
record layout instead of answering md5-only. The record is self-checking:
`norb` is written inside it and fixes the length of the `pnl` block, so the
layout must consume exactly the bytes gfortran's length marker claims —
`read_records` raises if it does not, which is what makes drift from the
Fortran loud rather than plausible. Integers and character fields are
treated as structure (no tolerance applies to `iz`, `nnlz`, `nlbcdu`); only
the `real*8` fields are compared numerically. 15 tests in
`CODE/tests/test_tape2n.py`; the full suite is 279 passed, 0 failed.

B2 against the baseline, with **zero structural differences** in every
record of every case — no integer, label, orbital count or mesh length
moved:

| case | records | reals | exact | max rel | rms |
|---|---|---|---|---|---|
| SnVIII | 2 | 41615 | 37926 | 2.41e-06 | 3.3e-08 |
| Forbid | 1 | 19907 | 19894 | 2.41e-06 | 3.4e-08 |
| Pr3 | 39 | 1066334 | 863016 | 4.01e-01 | 3.9e-04 |

`TAPE2N` moves by *less* than `OUT36` — 2.4e-06 against 1.2e-04 — and the
reason is worth recording: the alpha²-scaled quantities (`zeta`, the
mass-velocity and Darwin terms) are printed by RCN but are **not** carried
in the record; RCN2 recomputes them. What the record carries is the mesh
(unchanged to the bit), `ru`, `pnl` and the eigenvalues `ee`, and the
eigenvalues are what moved. Pr III's 4e-1 outlier is `pnl[21764]`, a
radial wavefunction in the far tail where the value itself is 6.6e-08; the
absolute move is 2.7e-08 against a function of order 1, so the relative
figure is meaningless there rather than alarming.

**Not done, deliberately:**

- **`CODE/RCN36K.exe` was not rebuilt in this section.** Replacing it
  promotes the new numbers, and per the hold point below that is the
  author's call. Done later, at the promotion step after Section 6, where
  all three executables were rebuilt and committed together.
- Section 9's documentation, which the plan schedules at the end. The
  `RCN36K.F` header entry (Vers. 36.4.2) was written now, because a change
  list that is written later is written from memory.

**Hold point (open question 3): promotion of the new `TAPE2N`.** The
numbers are above; whether they become the benchmark for RCN2 and RCG is
not this section's decision.

**Sequencing note for Sections 5 and 6 — read before starting them.**
`RINF` is shared: `RCN2K.F` uses it at 427 and 1078 and `rcg11kd.f` at
11726, 13665, 14581 and 14724. Its value has now been updated, so a
straight Stage A in those files cannot be byte-identical against a legacy
build — the §3 protocol assumed each file's constants were its own, and
`RINF`, `ALPH2H` and `PICON` are not. Only the eleven names RCN uses were
updated for exactly this reason; `ALPINL`, `ABOHR8`, `EV2CM`, `EV2CMK`,
`RYEV`, `RYKEV`, `GFE1`, `GACON`, `GACONX`, `PICONA` and `FOURO3` are still
at legacy values and are Sections 5 and 6 to update.

The mechanism for the shared three: every updated entry in `fundconst.f`
keeps its legacy value on the line above it as `C     LEGACY: PARAMETER
(...)`. A Section 5 or 6 Stage A is run from a build with those lines
uncommented and the CODATA ones commented out — a throwaway verification
build, never committed — which restores bit-exactness for the substitution
step. Then they go back and the value diff is read on its own.


### Section 5 — RCN2 (`RCN2K.F`) — **DONE 2026-08-20**

**Four** routines, not three: the 4/3 at line 618 is in `DIP`, which the §1
inventory did not list. 13 live sites — nine `109.73731D0` (`RINFKK`), two
`109737.31D0` (`RINF`), one `2.66246D-5` (`ALPH2H`), one `1.333333333D0`
(`FOURO3`). The pi at 2318 is inside commented-out code and is not a site, as
Section 3 already recorded.

`ZETA2`, `DIP`, `SLI2` and `G5INP` each carry `INCLUDE 'fundconst.f'` after
their `IMPLICIT` card and nothing else. Collision check re-run over
`RCN2K.F` for all 38 `fundconst.f` names: zero hits.

**Two stages, not three.** Every literal here already carried a `D0` suffix,
so the §3/§4 precision problem does not arise.

| Stage | What changed | Result |
|---|---|---|
| A | literal → parameter, `fundconst.f` at legacy values | **byte-identical** |
| B | `FOURO3` → `4.D0/3.D0`; `RINF`, `RINFKK`, `ALPH2H` already CODATA | values only |

The Stage A build used the §4 sequencing mechanism: `legacyize.py` in the
work area rewrites `fundconst.f` by activating each `C     LEGACY:` line and
commenting out the CODATA one below it — nine entries — for a throwaway build
that was never committed. Before any of it, the unmodified `RCN2K.F` was
rebuilt and shown to reproduce the Section 2 baseline: `ING11` md5-identical
on all three cases, `OUT2` identical apart from timing.

Stage A is byte-identical on `SnVIII`, `Pr3` and `Forbid`. `ING11`
md5-identical; `OUT2` identical apart from its timing lines.

**A third spelling of "timing line".** `OUT2`'s is a trailing ` SEC`, which is
neither RCG's `TIME=`/`FINISHED … AT`/`DIAGTIME` nor RCN's `t(run)=`/`t(job)=`,
so `cmp_num.py`'s default ignores do not cover it and it was passed as
`--ignore 'SEC *$'`. It was **not** added to `DEFAULT_IGNORES`, because 266 of
Pr III's 267 such lines are `I5,I5,I5,3F10.3,' SEC'` — ignoring the line would
throw away three structural counters to suppress three clock readings. Those
counters were checked separately, by diffing `OUT2` with only the time
*fields* blanked: **no counter moved in any case**. Whether the tool should
grow field-level ignores is left to the author; see §6, question 4.

**Stage B — the RCN2-only shift** (RCN held at legacy, i.e. both stages read
the Section 2 `TAPE2N`):

| case | file | numbers | max rel | rms | shape |
|---|---|---|---|---|---|
| SnVIII | ING11 | 83 | 5.18e-08 | 8.0e-09 | 81 exact |
| SnVIII | OUT2 | 98 | 6.78e-08 | 1.0e-08 | 95 exact |
| Forbid | ING11 | 56 | 5.18e-08 | 6.9e-09 | 55 exact |
| Forbid | OUT2 | 66 | 5.18e-08 | 6.4e-09 | 65 exact |
| Pr3 | ING11 | 4442 | 6.59e-06 | 1.0e-07 | 4397 exact |
| Pr3 | OUT2 | 21913 | 6.59e-06 | 6.3e-08 | 20883 exact |

The floor is **5.18e-08**, which is exactly the relative error of the legacy
`109737.31` / `109.73731` against CODATA 2022 — the whole of this section's
effect, and nothing else. `FOURO3` moves by 2.5e-10 and is invisible.

Everything above that floor is one quantum of a print field, checked rather
than assumed. All 45 Pr III `ING11` deviations resolve into three groups:
39 at abs 1.050e-01, which is 5.18e-08 × the 2.03e+06 kK `EAV` they sit on;
5 at abs 1.000e+01, in a column whose every entry ends in zero; and one at
abs 1.000e-04 on an `F9.4` field, a value sitting on a rounding boundary that
a 7.9e-07 shift tipped over. No structural field, label or line count moved
anywhere.

Note when reading `ING11` numbers out of `cmp_num.py`: the parameter cards are
`5(F9.4,I1)` and `7(F9.4,I1)`, a value with a one-digit precision code butted
straight onto it with no separator, so the tool reads value and code as one
number. The *relative* figures survive this — numerator and denominator both
scale by ten — but the absolute ones are ten times the true last-place unit,
which is why a single last-digit flip is reported as `abs 1.000e+01`.

**Cumulative** (RCN Stage B2 `TAPE2N` into RCN2 Stage B, against the Section 2
baseline), which is what RCG will actually receive:

| case | ING11: outside | max rel | OUT2: outside | max rel |
|---|---|---|---|---|
| SnVIII | 3 of 83 | 1.36e-06 | 3 of 98 | 1.36e-06 |
| Forbid | 1 of 56 | 9.95e-07 | 1 of 66 | 9.95e-07 |
| Pr3 | 50 of 4442 | 1.10e-04 | 1113 of 21913 | 6.60e-06 |

Pr III's 1.10e-04 is a single parameter printed as `9100` moving to `9101` —
one unit in the last of four significant digits. `ING11` carries parameters to
four significant figures, so that is the coarsest quantum in the file and the
largest move the format can express short of a real change. It is consistent
with RCN's own 1.2e-04 `OUT36` shift; nothing is amplified in RCN2.

**`ZETA2` was not tested**, because it is never used in practice: nothing
calls it except a `ZET` card in `IN2`, and its result is printed to `OUT2`
and goes nowhere else. Its two sites (`ALPH2H`, `RINF`) are textually exact
substitutions and the collision check clears them. The probe confirms the
rest: `RINFKK` and `FOURO3` are live in all three cases, `RINF` in Pr III
only, through the `SLI2` site at line 1078.

**Also worth knowing:** Norton deletes a freshly built `RCN2K` binary within
about a minute, under any filename, anywhere outside `CODE/` — an escalation
of the known false positive, which until now only ate *copies*. A binary
built into `CODE/` survives. Section 6 should build RCG's verification
executables there, under a temporary name, and delete them afterwards.

**Not done in this section, deliberately:** `CODE/RCN2K.exe` was not rebuilt
here, for the same reason as `CODE/RCN36K.exe` in Section 4 — see "Promotion"
after Section 6, where all three executables were rebuilt and committed
together.

### Section 6 — RCG (`rcg11kd.f`) — **DONE 2026-08-20**

**Five** routines, not six, and **26** sites, not 21: `LVDIST` has no literal
(the §1 count predates the §3 recount, which already said 26). `ENERGY`,
`SPECTR`, `WNDIST`, `BORN` and `RCOEFF` each carry `INCLUDE 'fundconst.f'`
after their `INCLUDE 'rcgpar.f'` — that include, not an `IMPLICIT` card, is
where RCG's `IMPLICIT REAL*8` comes from. Collision check re-run over
`rcg11kd.f` for all 38 names: zero hits.

Before anything, the unmodified `rcg11kd.f` was rebuilt into the work tree and
shown to reproduce the Section 2 baseline exactly on all three cases.

| Stage | What changed | Result |
|---|---|---|
| A | 26 literals → parameters, `fundconst.f` at legacy values | **byte-identical** |
| B | values → CODATA 2022 | values only |
| B3 | `RINFL0`→`RINF`, `ALPINL`→`ALPHIN`, `PICONA`→`PICON`, `GACONX`→`GACON`, now equal by definition | **byte-identical** |

Stage A and Stage B3 are byte-identical on all three cases: `OUTGINE`,
`TAPE2E` and `genealogy.txt` md5-identical, `OUTG11` identical apart from the
three timing lines. CLAUDE.md §1 applies to them in full.

**Stage B — the RCG-only shift**, with **zero structural differences**
anywhere: no line count, label or field width moved in any case.

| case | file | numbers | outside | max rel | rms |
|---|---|---|---|---|---|
| SnVIII | OUTG11 | 53148 | 744 | 3.390e-03 | 4.2e-05 |
| Forbid | OUTG11 | 4302 | 270 | 1.404e-03 | 1.6e-04 |
| Pr3 | OUTG11 | 9586062 | 525 | 4.545e-03 | 3.5e-06 |

`OUTGINE`, `TAPE2E` and `genealogy.txt` are **byte-identical in every case**.
This section's constants live entirely in the radiative and plotting code; the
levels, the eigenvectors and the RCE handoff do not see them.

Every deviation is one last-place unit of a print field, and the underlying
shifts are the constants: gf(E1) 1.87e-05 (`GFE1`), gA(E1) 8.6e-05 (that plus
`GACON`), M1/E1 5.8e-05 (two powers of `ALPINL`), E2/E1 6.7e-04 (two powers of
`ABOHR8`, whose 3.35e-04 was the largest single error in `fundconst.f`). The
Pr III summary block shows them unrounded: `SUMGF` 2957.9643 → 2958.0196 =
1.87e-05, `SUMGAR` 8.7966E+13 → 8.7974E+13 = 9.1e-05. The 1e-3 maxima are all
three-decimal `log gf` fields, which is §0 again: one unit there is 3.4e-03 in
gf and says nothing about the arithmetic.

**Cumulative** (Section 4 `TAPE2N` → Section 5 RCN2 → this RCG), against the
Section 2 baseline:

| case | OUTG11 outside / total | max rel | OUTGINE outside | max rel |
|---|---|---|---|---|
| SnVIII | 2526 / 53148 | 3.39e-03 | 181 / 542 | 1.54e-06 |
| Forbid | 271 / 4302 | 1.40e-03 | 0 / 99 | 0 |
| Pr3 | 167087 / 9586062 | 2.00e+00 | 472 / 15942 | 1.10e-04 |

`OUTGINE` moves only here and never in Stage B, so that shift is upstream and
none of it is this section's. Pr III's 2.00e+00 is not a finding: the extreme
relative deviations are eigenvector components printed `F8.5` whose values are
1e-5 — one or two units of the last decimal of a field whose own floor is
1e-5 — and every one has an absolute deviation of 5e-5 or less; with
`--atol 5e-5` the count falls from 167087 to 10621. `TAPE2E` is md5-identical
even cumulatively, and Pr III's line-selection counter (`2564994 LINES
OMITTED`) is unchanged in every run of every stage.

**Coverage — only 6 of the 26 sites are exercised.** Each site was perturbed
by 1% in turn, rebuilt and run; the silent ones were then re-probed against
all four output files rather than `OUTG11` alone. Live: `GFE1` (13187, all
three cases), `ALPHIN` (13188, Forbid, M1 only), `PICON` and `ABOHR8` (13190,
Forbid, E2 only), `GACON` (13502, SnVIII and Pr3; 13849, all three). The other
twenty — the `ENERGY` sites, the keV and plotting path in `SPECTR`, five more
`GA` branches, `WNDIST`, `BORN` and `RCOEFF` — are not reached by any deck
here. Those substitutions are textually exact and the collision check clears
them, which is the standing of Section 5's `ZETA2` sites; exercising them needs
a Born run, a plot run and a `TEXCIT` run, and no deck here has one.

**Defect 0a.1 is fixed but unexercised.** 13397, the `TEXCIT`-weighted
line-selection threshold, is dead in all three cases because none of them sets
`TEXCIT`. So the line-count change this section was told to expect did not
occur anywhere, and cannot be demonstrated with the present decks. The
reasoning is recorded at that line in `rcg11kd.f`.

**A caveat about `legacyize.py` from here on.** It now covers 17 entries, but a
legacyized build is no longer a reproduction of the *original* program: the
four sites collapsed in Stage B3 used coarser roundings (`109737.0`, `137.04`,
`3.14159`, `0.6698`) than the names they were folded into, and legacyize
restores the finer legacy value at all four. Everything else still round-trips.
For an original-output reproduction start from `rcg11kd.f.pre` in the work
tree, which is the file at commit 2ed5631.

**A `cmp_num.py` defect found and fixed on the way in.** Pr III's `OUTG11`
carries true-jj basis-state labels such as `232E1518`. The tokenizer read that
as 232e1518 = +Inf, two Infs gave Inf−Inf = NaN in `deviation`, and the whole
comparison died — the first time the tool had been pointed at a Pr III
`OUTG11`, since Sections 4 and 5 compare only `OUT36`, `OUT2` and `ING11`. An
exponent is now capped at three digits, which is all an `E` or `D` edit
descriptor ever writes, so `232E1518` tokenizes exactly as the `230D1518`
beside it always has. Regression test added; the suite is 280 passed.

**Not done in this section, deliberately:** `CODE/rcg11kd.exe` was not rebuilt
here, for the same reason as Sections 4 and 5 — see "Promotion" after this
section, where all three executables were rebuilt and committed together. Section 9's
documentation is still scheduled at the end; the `rcg11kd.f` header entry was
written now.

Work tree and full numbers: `F:\COWAN\fundconst_work\rcg\README.txt`.

Cases: `WORK/`, `Pr3/`, `Forbid`.

### Promotion of the rebuilt executables — **DONE 2026-08-20**

Sections 4, 5 and 6 each left `CODE/*.exe` alone on purpose, so that the
committed executables and the committed sources went out of step together
rather than one file at a time. Open question 3 settles it: the new `TAPE2N`
is the benchmark from here on, and the executables are rebuilt and committed
**together**, which is what keeps a mixed-alpha build from ever being
distributed — RCN's `zeta1` on the modern alpha while RCG's M1/E1 ratio still
carried `137.04` was the hazard, and it existed only between Section 4 and
Section 6.

`CODE/RCN36K.exe`, `CODE/RCN2K.exe` and `CODE/rcg11kd.exe` rebuilt from the
current `FOR/` sources with the Section 2 flag sets, straight into `CODE/`
(the Norton-safe location, per Section 5). `CODE/RCE20KD.exe` is untouched:
Section 8 has nothing in it.

Verified before committing, in two steps, in `F:\COWAN\fundconst_work\promote`
and `promote2`:

1. A build of the three programs into a scratch `bin/`, run end-to-end
   (RCN → RCN2 → RCG) on `SnVIII`, `Forbid` and `Pr3`, reproduces every
   recorded stage output: `TAPE2N` against Section 4's `stageB2`, `ING11`
   against Section 5's `rcn2/cum`, `OUTGINE`, `TAPE2E` and `genealogy.txt`
   against Section 6's `rcg/cum` — **md5-identical in all three cases** —
   and `OUT36` and `OUTG11` identical apart from their timing lines.
2. The executables actually committed were then run on all three cases and
   match that verified build, again apart from timing. A gfortran build is
   not byte-reproducible, so this is checked on output, not on the binary.

The one intended new difference anywhere is `OUT2`'s timing layout; see open
question 5. Test suite at the commit: 280 passed, 252 skipped.

### Section 7 — RCG conversion-coefficient routines (new code)

Only RCG needs these. Two new routines, per the author's specification:

- `GAFROMSCOEF(ITYPE)` — coefficient for gA from S; `ITYPE` = 1/2/3 for
  E1/M1/E2.
- `GFFROMGACOEF()` — gf from gA; type-independent.

Both `DOUBLE PRECISION FUNCTION`.

Four things to get right:

1. **`ITYPE`, not `TYPE`.** Under `IMPLICIT REAL*8 (A-H,O-Z)` a dummy named
   `TYPE` is a *real*. It would silently work for 1.0/2.0/3.0 and be a trap
   forever after. `TYPE` is also a keyword in later standards.
2. **Thread safety.** CLAUDE.md §2: `-fno-automatic` accompanies `-fopenmp`,
   so locals *and the function result variable* are static. These functions
   must not be called from inside a parallel region. Call them once during
   setup and cache the values in locals or `COMMON` — faster anyway, and it
   sidesteps the hazard entirely.
3. **Naming.** §6 asks new routines to group with their family.
   `GAFROMSCOEF` / `GFFROMGACOEF` are the author's names and will be used;
   noting only that the established pattern is a short family prefix (`CFC*`,
   `LSJJ*`).
4. **The coefficients themselves**, derived from
   `WORK/A-S-f-conversion_formulas.docx` and verified algebraically against
   the code's existing structure (§0b). With sigma in cm-1 and S in a.u.:

   - gf(E1) = `2/(3·AUCM)` · sigma · S = 3.0375568e-6 · sigma · S
   - gf(M1) = (alpha²/4) × the E1 coefficient
   - gf(E2) = (pi²·a0²·sigma²/5) × the E1 coefficient
   - gA (all types) = `8·pi²·alpha²·a0·c` · gf · sigma² = 0.66702518 · gf ·
     sigma²

   Note on the docx: the general Grant factor `Np` does not survive text
   extraction from the equation objects unambiguously. It was reconstructed as
   `Np = 2^(2p+2)·(p+1)·(2p+1)·pi^(2p+1) / (p·[(2p+1)!!]²)` and verified to
   reproduce the explicit E1, E2 **and** E3 formulas in the same document
   exactly. If a general-p routine is ever wanted, that is the form; for
   E1/M1/E2 the explicit closed forms above are used directly, so the
   reconstruction is not load-bearing here.

Acceptance here is different: this is new code with an independent oracle. The
routines must reproduce the Stage B constants to full double precision.
Unit-test them in `CODE/tests/`.

### Section 8 — RCE (`RCE20KD.F`)

No constants. Verify by inspection that it stays byte-identical, then close it
out. Included only for completeness of the top-to-bottom pass.

### Section 9 — Documentation

- ~~**The §0c `gA`-is-primary reminder, in all three required places**~~ —
  **done 2026-08-20**, ahead of the rest of the plan. See §0c.
- `FOR/LSJJ/README.txt` — long-form account, appended (§6: add to it, do not
  start a parallel document). Include the §0 quantization analysis and the
  seven-line verification table, so the investigation does not have to be
  repeated a third time.
- `rcg11kd.f` header change list — entry.
- `README.md` — last, after everything else is complete (§3).

---

## 5. Effort

| Section | Estimate |
|---|---|
| 1 comparator tool | 2–3 h — **done** |
| 2 baselines | 1 h (mostly Pr III runtime) — **done** |
| 3 `fundconst.f` | 1 h — **done** |
| 4 RCN | 3–4 h (the single-precision hazards) — **done** |
| 5 RCN2 | 1 h — **done** |
| 6 RCG constants | 1.5 h |
| 7 RCG coefficient routines | 2–3 h |
| 8 RCE | 15 min |
| 9 docs | 1.5 h (incl. §0c reminder x3) |

**~2 working days**, roughly 25% editing and 75% verification. Sections 1–3
are prerequisites; 4 through 8 are strictly sequential, because each consumes
the previous stage's output.

## 6. Open questions

1. ~~Widen the `LOG GF` print field?~~ **Resolved: no.** The field stays as it
   is; the reminder of §0c is the deliverable instead.
2. ~~**CODATA values** for `h`, `me`, `eps0`, `mu0`, `muB` are cited from
   memory in this document.~~ **Resolved 2026-08-20.** Every value in
   `FOR/fundconst.f` was transcribed from the NIST complete listing
   (`physics.nist.gov/cuu/Constants/Table/allascii.txt`, 2022 adjustment) and
   every derived quantity computed from those digits. Nothing in this document
   is load-bearing any more; the file is.
3. ~~**Promotion of the new RCN `TAPE2N`** to benchmark status.~~
   **Resolved (author, 2026-08-20): promote it.** The new `TAPE2N` is the
   downstream benchmark from here on, and `FOR/RCN36K.F` and
   `CODE/RCN36K.exe` are to be rebuilt and committed. One qualification,
   raised and accepted: between that commit and the end of Section 6 the
   shipped package is internally inconsistent in the fine-structure
   constant, because RCN's `zeta1` now uses the modern alpha while RCG's
   M1/E1 ratio still carries `137.04`. The fix is sequencing, not a
   different decision — the executables are rebuilt and committed together,
   once Section 6 lands, so no mixed-alpha build is ever distributed.
   **Carried out 2026-08-20**: `RCN36K.exe`, `RCN2K.exe` and `rcg11kd.exe`
   rebuilt and committed in one commit with the sources that produced them.
   See "Promotion" after Section 6.
4. ~~**Stage A of Sections 5 and 6** cannot be byte-exact against a legacy
   build for the constants RCN has already updated.~~ **Resolved (author,
   2026-08-20)**: the `C     LEGACY:` throwaway-build mechanism is accepted,
   and Section 5 used it successfully — `legacyize.py` in the work area, nine
   entries reverted, Stage A byte-identical on all three cases.
5. ~~**Should `cmp_num.py` grow field-level ignores?** New, from Section 5.
   `OUT2`'s timing lines are `I5,I5,I5,3F10.3,' SEC'`: three structural
   counters and three clock readings on one line, so a line-level ignore
   cannot suppress the clock without also discarding the counters.~~
   **Resolved (author, 2026-08-20): the comparator is not to change; RCN2's
   output is.** The clock readings now go on a line of their own, tagged
   `TIME=`, which `DEFAULT_IGNORES` already covers, and the three counters
   keep a line to themselves and are compared. Nothing reads `OUT2`, so the
   extra record is harmless. Effect on the Pr III `OUT2` comparison: 22708
   numbers compared with the default ignores, against 21913 under
   `--ignore 'SEC *$'` — the counters are back in the comparison and no
   by-hand check is needed. The change is deliberately unremarked in the
   source and in the documentation.
