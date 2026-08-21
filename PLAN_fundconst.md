# Plan — factor all fundamental constants into `FOR/fundconst.f`

Status: **complete, 2026-08-21.**  Written 2026-08-20.  Every numbered
section is done — 0c, and 1 through 13 — and the executables have been
promoted.  `CODE/RCN36K.exe`, `CODE/RCN2K.exe` and `CODE/rcg11kd.exe` are
built from the committed sources and reproduce the Section 13 reference on
all five regression cases; `CODE/RCE20KD.exe` is untouched, Section 8 having
nothing in it.  Record: `F:\COWAN\fundconst_work\promote3\README.txt`.

**Sections 10, 11 and 12 were added 2026-08-20, after Section 7, on the
author's instruction.**  Section 7 delivered `GAFROMSCOEF` and
`GFFROMGACOEF` and then left them uncalled, with `SPECTR` still folding the
same coefficients into `PSQ` and `GACON` by hand.  That is two independent
derivations of the same physics in one file, which is exactly what this plan
exists to remove.  A survey done at the same time found that Sections 4–6 also
*missed* three sites, one of them worth 5.2e-4 — larger than any constant
update this plan has made — and that the Python toolbox carries its own
hard-coded copy of the conversion coefficients.  The three new sections close
all of it.

**Section 13 was added 2026-08-20, on the author's instruction, and that
instruction is the standing rule of this plan:** *every fundamental constant
used anywhere in the package must come from `fundconst.f`.* Not most of them,
not the ones whose error is large enough to be worth a rebuild — every one.
A constant left inline is a constant that a future change of the values
cannot reach, and reaching all of them is the entire purpose of the file.
Section 10 quantified four that are still written inline and deferred the
decision to the author; the decision is made, and Section 13 folds them in.
Nothing in this document may say that a physical constant is to stay inline.

**Execution order, as carried out:** Sections 8 through 13, then Section 9,
then the promotion.  Section 9 kept its number and its place in this file so
that file order stays numeric order; only its execution was deferred.

Scope: replace every hard-coded physical constant in the Fortran with named
parameters from a single include file, update those values to CODATA 2022,
add exact E1/M1/E2 conversion coefficients to RCG, and make every consumer in
the package — Fortran and Python — obtain its coefficients from that one
place, so that nothing can escape a future change of the constants.

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

### Final promotion of the executables — **DONE 2026-08-21**

The three executables committed after Section 6 were five sections out of
date: Sections 7, 10, 11 and 13 all moved `rcg11kd.f`, and Section 10 moved
`RCN2K.F` as well. Work area `F:\COWAN\fundconst_work\promote3\`, full
record in its `README.txt`.

**One thing had to be undone first.** `CODE/rcg11kd.exe` in the working tree
was a byte-for-byte copy of `sweep/binS11/rcg11kd.exe` — the Section 11
build, copied in during that section and stale by the time Section 13 landed.
It was reverted before Sections 10–13 were committed, so that no half-promoted
binary reached the history at all.

Built with the Section 2 flag sets, gfortran 13.2.0 (MinGW-W64, Strawberry);
the only warning anywhere was the expected `-fno-automatic overwrites
-frecursive`. `RCN2K.exe` was built straight into `CODE/` under a temporary
name and renamed there, because Norton deletes a freshly built `RCN2K` binary
from anywhere under `fundconst_work` — the Section 5 and 10 workaround.
`CODE/` holds no scratch binary now; checked.

Verified in three steps:

1. The scratch build run **end to end — RCN → RCN2 → RCG** — on all five
   cases, against `gfac/stageB2`: **0 of 40 file comparisons differ**
   (`TAPE2E`, `OUTGINE`, `ING11`, `TAPE2N`, `genealogy.txt` raw; `OUT36`,
   `OUT2`, `OUTG11` less timing lines). This is stronger than Section 13's
   own check, which seeded RCN and RCN2 output and ran only RCG on the
   strength of the object-file identity proof: here they were actually
   re-run, and `TAPE2N`, `OUT36`, `OUT2` and `ING11` come back identical
   anyway.
2. The executables **actually placed in `CODE/`**, run on the same five
   cases: 0 differ against the scratch build and 0 against `gfac/stageB2`.
   A gfortran build is not byte-reproducible across invocations, so this is
   checked on output rather than on the binary, as the Section 6 promotion
   was.
3. **The CFP decks rebuilt by the promoted `rcg11kd.exe`**, following
   `make_cfp.bat` by hand: `FOR072`, `FOR073` and `FOR074` all come back
   md5-identical to the decks in use. Not asked for by this plan and worth
   having — the decks are unformatted and machine-built, and a user who
   recompiles RCG is told to rebuild them. This recompilation does not
   require it, and now that is established rather than assumed.

`CODE/RCE20KD.exe` is untouched. The plan is complete.

### Section 7 — RCG conversion-coefficient routines (new code) — **DONE 2026-08-20**

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

**Done.** Both functions are at the end of `FOR/rcg11kd.f`, built entirely
from `fundconst.f`, with no `COMMON`, no unit, no `SAVE`, no `DATA` and no
state of any kind — which is what lets them be compiled and tested on their
own. The usage is

```
gA = GAFROMSCOEF(ITYPE) * S * SIGMA**(2*p+1)
gf = GFFROMGACOEF()     * gA / SIGMA**2
```

sigma in cm-1, S in a.u., gA = g(upper)·A in s-1, gf = g(lower)·f. The power
of sigma is 2p+1 — 3 for E1 and M1, **5 for E2**. That is not an
inconsistency in the coefficient; it is what makes the coefficient a
constant.

| | | value returned | handbook |
|---|---|---|---|
| `GAFROMSCOEF(1)` | E1 | 2.02612688449478421D-06 | 2.0261e18 |
| `GAFROMSCOEF(2)` | M1 | 2.69735002205593345D-11 | 2.6974e13 |
| `GAFROMSCOEF(3)` | E2 | 1.11995002964536993D-22 | 1.11995e18 |
| `GFFROMGACOEF()` | all | 1.49919378619413113D+00 | 1/GACON |
| any other `ITYPE` | — | 0 | — |

**The oracle held.** Martin & Wiese publish gA(E1) = 2.0261e18·S/λ³,
gA(M1) = 2.6974e13·S/λ³ and gA(E2) = 1.11995e18·S/λ⁵ with λ in Å; since
sigma = 1e8/λ those are the three coefficients times 1e24, 1e24 and 1e40, and
they agree **to every digit the handbook prints**. This is the only
independent check on `fundconst.f` anywhere in the package — alpha, a0, c and
the hartree all enter, and getting any of them wrong moves one of the three.

The compiled functions were checked against the same expressions evaluated in
Python from the CODATA 2022 digits and agree **bit for bit**, so the tests
assert float equality rather than a tolerance: nothing rounds anywhere between
the CODATA digits and the answer.

The four things the section asked to get right: `ITYPE` is declared
`INTEGER*4` explicitly and a test asserts both the name and the declaration,
since neither shows in the returned values; the non-re-entrancy under
`-fno-automatic` is stated in the routine header, the file header and the test
file, and nothing calls either function yet so nothing can be violating it;
the author's names were used as given, noting only that they carry no family
prefix; and a test enforces that the only literals in the arithmetic are the
exact small integers 0, 1, 2 and 5.

`CODE/tests/test_coefficients.py`, nine tests in two layers — a Python layer
that carries the handbook comparison and always runs, and a layer that
compiles the two routines out of `rcg11kd.f` and skips without gfortran.
Suite: **289 passed, 252 skipped** (280 before).

**Not done in this section.** Nothing calls them: `SPECTR` still folds the
same coefficients into `PSQ` and `GACON` as it always has. Rewiring it changes
a computed result and needs its own §1 argument; this section did not ask for
it. Because the addition is dead code it cannot move any existing output, and
no case was re-run. **That is now Section 11**, added on the author's
instruction — leaving the two derivations side by side is precisely the
duplication this plan exists to remove, and “dead code” is the wrong resting
place for the only routines in the package that have an independent oracle. `CODE/rcg11kd.exe` was **not** rebuilt, for the same
reason Sections 4–6 left the executables alone — **a final promotion is owed**
once Sections 8 and 9 close out.

Record: `F:\COWAN\fundconst_work\coef\README.txt`.

### Section 8 — RCE (`RCE20KD.F`) — **DONE 2026-08-20**

No constants. Verify by inspection that it stays byte-identical, then close it
out. Included only for completeness of the top-to-bottom pass.

**Done, and the verification was made mechanical rather than by eye.** The
Section 10 scanner was written early and run first over the whole RCE family
— `RCE20KD.F`, `RCE20KD_PAR.F`, `RCEJJ2E.F`, `LSJJCOM.F`, `LSJJID.F`,
`LSJJPAR.F` — flagging every non-comment numeric literal outside `FORMAT`
statements and outside character strings. Result: **not one real literal with
three or more significant digits anywhere in the six files.** The most precise
reals in the whole of RCE are `0.99D0` (a step-damping clamp in the
least-squares iteration), `0.001E0` (a divide-by-zero guard in the same
place), `0.5D0`, `2.D0` and `60.D0` (seconds to minutes). All of the many
three-and-four-digit *integer* hits are statement labels, `FORMAT` labels,
array dimensions, the `MAGJJ=20260802` deck stamp, and the eigenvector-distance
thresholds `1000000` / `58000000`. None is a physical constant.

Two small corrections to what §1 asserted, both in RCE's favour:

- §1 said "the only decimals are in comments". Not literally true — there are
  the six algorithmic reals above — but the *conclusion* is right: none of them
  is a fundamental constant, and none belongs in `fundconst.f`. The decimals
  that do appear in RCE's comments (`191.931`, `1.573E+51`, `0.95`, `0.76`,
  `3.9 S`, `0.69 S`) are fit diagnostics and threading timings, not physics.
- §1 credited RCE with 0 sites in 178 routines package-wide; the scan confirms
  0 and extends the same result to `RCEJJ2E.F` and the three `LSJJ*.F`
  includes, which §1 listed separately but by the same argument.

`RCE20KD.F` and its includes are untouched in the working tree (`git status`
clean for all six), so byte-identity of the source is exact rather than
inferred, and `CODE/RCE20KD.exe` needs no rebuild for this plan. RCE therefore
gets **no `INCLUDE 'fundconst.f'`** — adding one would be a compile-time
dependency purchased for nothing.

**One consequence to state plainly, because it is easy to mistake for a
regression.** RCE's numbers still move, even though RCE did not change: it
reads `TAPE2E` from RCG, and Section 6 Stage B changed the constants that
produce `TAPE2E`. Any RCE output that shifts after this plan lands is
inherited from upstream, not introduced here. That is the expected behaviour
of a pure least-squares fitter, and it is the reason RCE needs no section of
its own beyond this one.

Record: `F:\COWAN\fundconst_work\rce\README.txt`.

### Section 9 — Documentation — **DONE 2026-08-21**

**Runs last**, after Sections 8, 10, 11 and 12; it keeps this number and this
position so that the file reads in numeric order. Add to its list, from the
sections added on 2026-08-20: the three missed sites and the autoionization
rate change of Section 10 (which is a *visible* change to `WORK/AI` output and
must be called out as such, not buried in a constants table); the removal of
the duplicate coefficient derivation in Section 11; and the Fortran-to-Python
tie of Section 12. The `rcg11kd.f` header change list gets an entry for each.

Three things Section 10 turned up that the documentation must carry, and
which are easy to lose because they are not value substitutions:

- The autoionization rates rise by **+4.81e-04**, not the +5.2e-04 predicted
  here, because `RINF` in the same expression pulls the other way. Quote the
  number that the runs show, not the one AURATE alone implies.
- **Four physical constants that Section 10 found still written inline** —
  the electron g-factor as `0.50116D0`, `1.50116D0` and `1.00232D0`, and the
  Maxwellian rate coefficient `2.1717D-8`. Section 13 folds them into
  `fundconst.f`; document the resulting shifts (7e-07 to 1.5e-05, and the
  printed Landé g-factors among them) with the rest of the Stage B tables.
  Once Section 13 has landed, the package documentation may state without
  qualification that every fundamental constant comes from `fundconst.f` —
  and it must not state it before, so write this section after Section 13,
  not alongside it.
- `legacyize.py` can no longer restore the original text of every site.

- ~~**The §0c `gA`-is-primary reminder, in all three required places**~~ —
  **done 2026-08-20**, ahead of the rest of the plan. See §0c.
- `FOR/LSJJ/README.txt` — long-form account, appended (§6: add to it, do not
  start a parallel document). Include the §0 quantization analysis and the
  seven-line verification table, so the investigation does not have to be
  repeated a third time.
- `rcg11kd.f` header change list — entry.
- `README.md` — last, after everything else is complete (§3).

**What was done.** Three places, in the order the section asked for, plus two
entries in the distribution's own file list.

- **`rcg11kd.f` header change list.** The fundconst entry already carried
  Sections 6, 7, 11 and 13; what it was missing was Section 10. Added: the
  sweep itself and why the §1 inventory could not have found these sites, the
  `2066` autoionization change written out in full with the before-and-after
  line and the 4.81e-04 arithmetic, `CSIGMA`/`PICSCN` recorded as
  **unexercised**, and the `RINF`-has-three-legacy-roundings caveat with what
  a bit-exact reversion now has to do instead of running `legacyize.py`. Two
  corrections to what the header already said: the site count was `34`, which
  Sections 10, 11 and 13 had made stale — it is 36 — and the claim that the
  g-factor is "the one place in this whole revision where a printed number the
  author reads directly changes" was true when it was written and stopped being
  true when Section 10 landed. It now says *second of the two*, with the
  autoionization rates named as the first. Comment-only; `-fsyntax-only` passes
  and every added line begins with `C` in column 1 and stops before column 72.

- **`FOR/LSJJ/README.txt`** — 442 lines appended under EVERY FUNDAMENTAL
  CONSTANT IN ONE FILE, in the file's own house style (§6: add to it, do not
  start a parallel document). It carries the §0 quantization analysis and the
  seven-line verification table in full, so the gA/log gf investigation does
  not have to be done a third time; both §0a defects; the design of
  `fundconst.f` and the three properties of it that are load-bearing; the
  two-stage protocol and why Stage A is worth more than it looks; the
  program-by-program account including RCE's *nothing*, and why RCE's numbers
  move anyway; the coefficient routines as the only independent oracle in the
  package; the sweep and the `2066` error; Section 13 and the rule; the Python
  side with both of its defects; and a closing table of every shift a reader
  comparing old output will meet.

- **`README.md`** — a new `### THE FUNDAMENTAL CONSTANTS` section between the
  additions to the manual and the hardware requirements, opening with *read
  this before comparing against an earlier release*: what changed and what
  did not, the shift table, the autoionization change in plain words, and a
  pointer at `FOR/LSJJ/README.txt` for the rest. Five entries added to
  SUMMARY OF THE MAJOR BUGS FIXED — the `2066` truncation, the `0.6698`
  dropped digit, RCN's two fine-structure constants, the inconsistent
  `gf(E2)`, and the factor-of-1000 `EV_PER_CM`. The version dates of RCN,
  RCN2 and RCG moved to Aug 2026; RCE's stayed at Jan 2026, because RCE did
  not change.

- **The addendum file list** gained `FOR/fundconst.f` and
  `CODE/toolbox/fundconst.py`, and the `FOR/LSJJ/` entry now says its
  `README.txt` is the long account of the constants work as well as of the
  coupling work.

Two things were deliberately *not* written into `README.md`: a pointer to
`PLAN_fundconst.md`, which is not in the distribution, and any figure for a
site the runs do not exercise presented as though it had been measured. The
unexercised sites are named as unexercised in `FOR/LSJJ/README.txt` and
counted in the `README.md` pointer sentence.

Suite unchanged by any of it: **309 passed, 252 skipped**.


### Section 10 — Sweep: the sites Sections 4–6 missed — **DONE 2026-08-20**

Sections 4, 5 and 6 worked from the §1 inventory, and the inventory was built
by searching for the *known* spellings of the constants. Anything spelled
unexpectedly was invisible to it. A re-scan — every non-comment numeric
literal with three or more significant digits in all four programs, read one
by one — turns up **three sites carrying physics that no earlier section
touched**, and one of them is worth more than any value change this plan has
made so far.

| Site | Literal | What it is | CODATA 2022 | Legacy error |
|---|---|---|---|---|
| `FOR/rcg11kd.f` `ENERGY`, autoionization rate | `2066.D0` | Ry/hbar in units of 1e13 s-1 | `2067.0686667586` | **−5.17e-04** |
| same line | `109735.3D0` | Rinf in cm-1 | `109737.31568157` | −1.84e-05 |
| same line | `6.28318D0` | 2·pi | `6.28318530717959` | −8.45e-07 |
| `FOR/rcg11kd.f` `SPECTR`, `CSIGMA` | `8.06662D-18` | 4·pi²·alpha·a0² in cm², the photoionization cross-section coefficient | `8.0672837256768D-18` | −8.23e-05 |
| `FOR/RCN2K.F` `G5INP` | `109.737D0` | Rinf in 1000 cm-1 | `109.73731568157` | −2.88e-06 |

The line in `ENERGY` is

```
AA(LP,I)=((SUM*UENRGY/109735.3D0)**2)*6.28318D0*2066.D0
```

i.e. `2·pi·|V|²` in rydbergs converted to 1e13 s-1 — Fermi's golden rule.
`2066` is a truncation, not a rounding, of `2067.07`; it is the single
largest constant error left anywhere in the package, and unlike defect 0a.1
it is **exercised**, by `WORK/AI` and `WORK/AI/DIEL`. Autoionization rates
printed by RCG were 0.05 % low.

`CSIGMA` is `4·pi²·alpha·a0²`, verified algebraically: it evaluates to
8.0672837e-18 cm² against the 8.06662e-18 in the source. It gates the
photoionization cross-section path (`SIGMA`).

The five sites are now

```
AA(LP,I)=((SUM*UENRGY/RINF)**2)*TWOPIC*AURATE
IF (G.EQ.GOLD.AND.H.EQ.HOLD) CSIGMA=PICSCN
EMAX=EAV(NC)-RINFKK*E(NCSPVS)
```

**Two departures from the specification above, both forced and both recorded
in `fundconst.f` at the point of change.**

- The name is `TWOPIC`, not `TWOPI`. The collision check this section asked
  for is the thing that found it: `RCN36K.F` subroutine `phshift` has a local
  `twopi`, and that routine includes `fundconst.f`, so a `PARAMETER` spelled
  `TWOPI` is an assignment to a parameter. `AURATE`, `PICSCN` and `HBARJ`
  were clean. This is the second time the collision check has changed a name
  (`ALPHFS`, `PICON` were the first), and it is the reason the rule at the
  head of `fundconst.f` says to re-run it.
- `HBARJ` is in **Section B**, not Section A. It is `PLANCK/(2.D0*PICON)`,
  `PICON` lives in Section B, and a `PARAMETER` may only use names defined
  above it. Section A would have needed a second spelling of pi in the one
  file that exists to prevent exactly that.

**The value change is +4.81e-04, not the +5.2e-04 this section predicted**,
and the difference is arithmetic rather than noise: `RINF` appears squared in
the denominator of the same expression, so its own −1.84e-05 pulls the other
way.

```
(2067.0686667586/2066) * (109735.3/109737.31568157)**2
                       * (6.28318530717959/6.28318)  =  1.000481354
```

which is what the runs show, off `SUMGAA`: `1.96195E+06 → 1.96290E+06`,
`6.48954E+06 → 6.49266E+06`, `9.05517E+06 → 9.05952E+06` in `AIDIEL` and
`1.58222E+15 → 1.58298E+15` in `AI`. `AIDIEL`'s `GM*FRBAR` and its whole
`T(KEV)` dielectronic-recombination table move with them,
`9.0552E+06 → 9.0595E+06`.

**Five regression cases, not three.** The two `WORK/AI` decks had never been
used as regression cases and were baselined first; without them this
section's largest change would have been unexercised. Both do reach the
`ENERGY` line — the plain `AI` deck prints `AA(UNITS OF 1/SEC)` tables
although it prints no `GAATOT`, which belongs to the DR summary that only
`AIDIEL` reaches.

The staging was the three-step form Sections 4–6 used:

| | | |
|---|---|---|
| `base` | committed source, CODATA values, legacy literals | — |
| `stageA` | sites substituted, parameters at LEGACY values | **byte-identical to `base`** |
| `stageB1` | values updated, temporary names still in place | the value change |
| `stageB2` | temporary names folded away = committed source | **byte-identical to `stageB1`** |

Stage A needed two names that are not in the committed file — `RINFL3 =
109735.3D0` and `RINFK7 = 109.737D0` — because `RINF` and `RINFKK` each now
stand for more than one legacy rounding and their single `LEGACY` line can
restore only one. Giving those two names the CODATA values in B1 makes
B1 → B2 a pure rename, and it is byte-identical, which is what proves the
fold. `base` → `stageB2` differs in **two files only**: `AI/OUTG11` and
`AIDIEL/OUTG11`. `TAPE2E`, `OUTGINE`, `ING11`, `TAPE2N` and `genealogy.txt`
are byte-identical in every comparison including that one.

`legacyize.py` alone can no longer restore the package's original text —
`RINF` now has three legacy roundings behind it and the script knows one.
The caveat recorded at the end of Section 6 is now one degree worse, and it
is written into `fundconst.f` under `RINF` rather than left in this plan.

**Coverage, stated rather than glossed.** `PICSCN` is unexercised: `CSIGMA`
gates the photoionization path and no deck in `WORK/` turns it on, so the
value was checked algebraically and by nothing else. `RINFKK` at `G5INP`
1718 is substituted and proved bit-exact by Stage A, but its value change
reaches no printed output in any of the five cases — `ING11` and `OUT2` are
byte-identical throughout.

**What the sweep found that this section did not list.** Four more literals
that *are* physical constants and were still written inline. They were **not**
changed *here*: each moves printed output, and widening a section that was
already carrying the autoionization change was the author's call to make.
It was made — see the standing rule at the head of this plan — and these four
are the subject of **Section 13**.

| Site | Literal | What it is | CODATA 2022 | Error |
|---|---|---|---|---|
| `rcg11kd.f` 1328 1329 8979 8980 9102 9103 15217 | `0.50116D0`, `1.50116D0` | the electron g-factor as (g−1)/2 and 1+(g−1)/2, in the Landé g-factor | `0.501159652180`, `1.501159652180` | +6.9e-07, +2.3e-07 |
| `rcg11kd.f` 10454 | `1.00232D0` | g−1, in the magnetic-moment matrix element | `1.002319304361` | +6.9e-07 |
| `rcg11kd.f` 14917 | `2.1717D-8` | 2·√pi·alpha·c·a0², the Maxwellian rate coefficient per unit collision strength, cm³/s, T in rydbergs | `2.1716669877D-08` | +1.5e-05 |

All four are at least an order of magnitude below the 5.2e-04 that made the
`2066` site worth finding, which is why deferring them by one section was
defensible; leaving them permanently would not have been. They are listed in
the new test in a group of their own, `PHYSICS_STILL_INLINE`, which is a
holding pen with an expiry date, not a category: **Section 13 empties it and
deletes it**, and the test's own docstring says so.

**The hole is closed.** `CODE/tests/test_literals.py` reads every non-comment
real literal of three or more significant digits out of `FOR/*.f` and
`FOR/LSJJ/` and fails on any that is not in an inventory of 198 tokens, each
with a reason: CODATA transcriptions, pure mathematics, Thomas–Fermi model
constants, quadrature weights, iteration control, print-field sentinels, sort
keys, one empirical fit, three `DATA` tables, and the four above. Three
tests — nothing unaccounted for; nothing *stale* in the inventory, so a
removed literal cannot be silently re-permitted; and the five Section 10
literals specifically have not come back.

Integers are out of scope on purpose, and the test says so at the top: in
this package an integer literal is a statement label, an array bound, a loop
limit or a `FORMAT` width — hundreds of them, not one of them physics — and
every physical constant here is written with a decimal point or an exponent,
`2066.D0` included. Pure mathematical constants are listed but are *not*
candidates for `fundconst.f`; they cannot go stale.

Full suite after the change: **292 passed, 252 skipped**.

Record: `F:\COWAN\fundconst_work\sweep\README.txt`.

### Section 11 — RCG calls the coefficient routines — **DONE 2026-08-20**

The point of Section 7. `SPECTR` currently derives, by hand,

```
PSQ=GFE1*(PMAX**2)
IF (II.EQ.0) PSQ=PSQ*(0.5D0/ALPHIN)**2
IF (II.EQ.2) PSQ=PSQ*((PICON*ABOHR8)**2)/5.D0
```

and then eight separate `GA=GACON*...` sites. `GAFROMSCOEF` derives the same
three type factors a second time, from the same constants, a few thousand
lines away. Both are correct today. Two derivations that must agree and are
not mechanically tied together will not stay agreeing.

**One new routine is needed.** `SPECTR`'s `PSQ` is a *gf-from-S*
coefficient, and Section 7 built only gf-from-gA and gA-from-S. Add

```
DOUBLE PRECISION FUNCTION GFFROMSCOEF(ITYPE)
      GFFROMSCOEF = GAFROMSCOEF(ITYPE)*GFFROMGACOEF()
```

with the same conventions, the same `ITYPE`-not-`TYPE` rule, the same
returns-0-on-a-bad-type rule, and the same thread-safety note. `gf =
GFFROMSCOEF(ITYPE) * S * SIGMA**(2p-1)` — sigma^1 for E1 and M1, sigma^3 for
E2. Extend `CODE/tests/test_coefficients.py` to cover it against the
handbook's `gf(E1) = 303.75568885954 S/lambda`, which is `GFE1`·1e8 and is
already the number `CODE/toolbox/radiate.py` carries.

Then, in `SPECTR`, once at set-up and outside every parallel region (§2 of
CLAUDE.md — these functions are not re-entrant under `-fno-automatic`):

```
      GFCF=GFFROMSCOEF(ITYPE)
      GACF=1.D0/GFFROMGACOEF()
      PSQ=GFCF*(PMAX**2)
```

`ITYPE` from `II` (`II=0` M1, `II=1` E1, `II=2` E2 — confirm that mapping at
the source before writing it down), and `GACF` replacing `GACON` at all eight
`GA=` sites. Caching in a local is what point 2 of the Section 7
specification already required, and it is faster than the parameter reference
was.

**Two things this section must get right.**

1. **The E2 scaling.** `SPECTR` writes the E2 factor as
   `(PICON*ABOHR8)**2/5` with a0 in units of 1e-8 cm, and then divides by
   `ERAS2**2` where `ERAS2 = 1.0D8/ERAS1` is the wavelength in ångström. The
   1e16 introduced and the 1e16 removed cancel *algebraically* but not *in
   floating point*: neither is a power of two. `GAFROMSCOEF(3)` uses a0 in
   cm and has no such pair. The header of `GAFROMSCOEF` already records that
   the two agree to one ulp of the intermediate. Rewiring therefore moves
   the last bit or two of the E2 numbers, and the same is true to a lesser
   degree of E1 and M1, where `GACON*GFE1*(1/GACON)` is not bit-for-bit
   `GFE1`.
2. **Acceptance, consequently, is not §1's bit-exactness of the arithmetic —
   but it is still §1's byte-identity of the output.** Nothing here changes a
   value; the deviation is a few units in the 16th significant figure, and
   the widest print field in `OUTG11` is nine. The standard for this section
   is: `TAPE2E`, `OUTGINE` and `genealogy.txt` md5-identical and `OUTG11`
   identical apart from timing, on `SnVIII`, `Forbid` and `Pr3` — `Forbid`
   being the only deck that exercises M1 and E2 at all, and therefore the
   one that matters. If a single field does move, it will be a knife-edge
   rounding in the last printed decimal; that is acceptable **only** if it is
   shown to be exactly that, one field at a time, with `cmp_num.py` at
   `--rtol 1e-12` clean everywhere else. Do not accept a diff you have not
   individually explained.

Consider reversing the E2 arrangement rather than preserving it, if that is
what makes the section clean: computing sigma in cm-1 throughout and dropping
the 1e8/1e16 pair entirely is closer to the physics and removes the
discrepancy at its source. It is also a larger change to a routine nobody
wants to disturb. Decide it explicitly, in the section, not while editing.

**Then delete the duplication.** Once `SPECTR` is on the routines, `GFE1`,
`ABOHR8` and the `(0.5D0/ALPHIN)**2` and `/5.D0` type factors have exactly
one consumer each — `GAFROMSCOEF`. `GFE1` and `GACON` stay in `fundconst.f`,
since the functions are built from them; `ABOHR8` should go, its one site
gone and `ABOHR` doing the work. Leaving an unused rounding-specific name in
`fundconst.f` invites its reuse.

Estimate: 3–4 h, most of it verification on `Forbid`.

**Decisions taken while executing the section.**

1. **The E2 arrangement is preserved, not reversed.** `SPECTR` keeps its
   ångström scaling and `GFFROMSCOEF(3)` is multiplied by `1.D16` at the one
   point `PSQ` is formed. Reversing it is not the local change the section
   above imagined: `PSQ` is consumed by seven `GF=` sites, which do divide by
   `LAMBDA**2`, but also by the `TEXCIT` line-selection weight, by both
   `BRNCH` sites, and by the unit-20 record that `WNDIST` reads back — none
   of which divides by anything. Putting sigma in cm-1 throughout would move
   the branching ratios and the total decay rates by 1e16 and then require
   each of those four paths to be corrected by hand. That is the total
   radiative decay rate calculation, which §3 of `CLAUDE.md` says stays as it
   is. The `1.D16` is commented at the point of use as the unit conversion it
   is.
2. **The `GA=` sites keep `GACON`; `GACF` was not introduced.** The section
   above asked for `GACF=1.D0/GFFROMGACOEF()` in place of `GACON` at the
   `GA=` sites. That is not a de-duplication: `GACON` at those sites is
   already a direct `PARAMETER` reference to the single definition in
   `fundconst.f`, not a second derivation of anything, and the substitution
   would replace it with a runtime reciprocal round-trip. Nor is it faster —
   a `PARAMETER` costs nothing at run time, so the claim above that caching
   "is faster than the parameter reference was" is simply wrong. It is
   bit-exact today only because `1/(1/GACON)` happens to return `GACON` for
   this particular value, which is luck, not a property. The duplication
   Section 7 exists to remove was the *type factors*, and `PSQ` was the only
   place they were written twice. (The section also says "eight `GA=` sites".
   There are eight, but only seven are in `SPECTR`; the eighth is in `BORN`.)
3. **`303.75568885954` is not `GFE1`·1e8.** `2e8/(3·219474.63136314)` is
   `303.75568352754544`, lower by 1.8e-8 relative — the published figure was
   computed from an older hartree. `CODE/perl_legacy/conv_out.pl` and
   `CODE/toolbox/radiate.py` both carry the published one. The new test
   therefore checks it to `rel=5e-8` and asserts explicitly that it is *not*
   an equality, with the reason written down beside the number. **Section 12
   inherits this**: rebuilding `radiate.py` on `fundconst.py` moves its
   `303.75568885954` to `303.75568352754544`, which is a real (if tiny) value
   change to that tool's output and must be called out as one.

**What was done.**

- `GFFROMSCOEF(ITYPE)` added at the end of `FOR/rcg11kd.f`, written as
  `GAFROMSCOEF(ITYPE)*GFFROMGACOEF()` so that the type factors are derived
  once in the file. It takes no constants of its own and includes nothing.
  Returns `3.03755683527545D-06`, `4.04385039225686D-11`,
  `1.67902212529227D-22`, and 0 on any other `ITYPE`.
- `SPECTR`'s inline `PSQ` derivation replaced by `GFFROMSCOEF`. The `II` →
  `ITYPE` map is `II=0`→2 (M1), `II=1`→1 (E1), `II=2`→3 (E2), confirmed at
  the source. The `IGEN.GE.5` branch — autoionization, which takes no E2 type
  factor and a factor 3 instead — the section above did not mention; it stays
  on `ITYPE=1` even when `II` is 2.
- `ABOHR8` removed from `fundconst.f`, its only site gone. `ALPHIN` also has
  no consumer left, but it is a plain CODATA-derived name rather than a
  unit-scaled variant, so it stays.
- `CODE/tests/test_coefficients.py` extended: six new assertions over
  `GFFROMSCOEF`, including that what it returns is bit-for-bit the product of
  the other two.

**Result.** The reassociation costs 0 ulp on E1, 2 on M1 and 4 on E2 —
measured, not estimated. All five sweep cases (`SnVIII`, `Forbid`, `Pr3`,
`AI`, `AIDIEL`) are **byte-identical** to the Section 10 baseline: `TAPE2E`,
`OUTGINE`, `ING11`, `TAPE2N` and `genealogy.txt` raw, `OUT36`, `OUT2` and
`OUTG11` less timing lines. No field moved, so no field needed explaining.
`Forbid` is the deck that carries M1 and E2.

### Section 12 — The Python side: one source for the toolbox constants — **DONE 2026-08-21**

`CODE/toolbox/radiate.py` carries its own table:

```
_GF = {E1: (1, 303.75568885954), M1: (1, 0.0040438504), E2: (3, 167.84224),
       M2: (3, 0.002235255),     E3: (5, 47.140897),    M3: (5, 0.000627579)}
_S  = {E1: (3, 2.0261269e18),    M1: (3, 269735e8),     E2: (5, 1.11995e18),
       M2: (5, 14909714e6),      E3: (7, 3.14441e17),   M3: (7, 4.1861e12)}
EV_PER_CM = 8.06554393734921
```

transcribed from `conv_out.pl:1650-1710`. These are hard-coded values of
exactly the kind this plan exists to remove: they are the *same physics* as
`GAFROMSCOEF`, they are consumed by `conv_out.py` and `calc_TP_M1E2.py`, and
a change to `fundconst.f` does not reach them. `radiate.py` is already the one
copy on the Python side — Section 12 makes it a *derived* copy instead of an
independent one.

Two candidate mechanisms; the second is the recommendation.

- **Parse `FOR/fundconst.f` at import.** Truest to "one source", and the
  worst engineering: `CODE/` is what the author copies to `C:\cowan\code`
  and puts on the path, `FOR/` is not distributed with it, and a Python
  module that fails to import when the source tree is absent is a bad trade
  for a value that changes once a decade.
- **A generated `CODE/toolbox/fundconst.py`,** carrying the CODATA 2022
  digits and the same derivations, with a test that holds it to the Fortran.
  This is the pattern Section 7 already established and validated:
  `test_coefficients.py` compiles the two routines out of `rcg11kd.f`,
  evaluates the same expressions in Python, and asserts **float equality**.
  Extend exactly that mechanism — the test compiles `fundconst.f` into a
  tiny driver that prints every Section A and Section B parameter in
  `E24.17` and asserts each equals its Python counterpart bit for bit. The
  test skips without gfortran, as the existing one does, so the toolbox
  stays importable anywhere while the tree stays provably consistent
  wherever it can be checked.

Then rewrite `radiate.py` on top of it:

- E1, M1 and E2 come from the same three expressions the Fortran uses. They
  must reproduce `GAFROMSCOEF(1..3)` to full double precision, scaled by the
  1e24, 1e24 and 1e40 that Section 7 already worked out — a test asserts it,
  and that test is the tie between the two languages.
- M2, E3 and M3 have no Fortran counterpart. Compute them from the general
  Grant factor reconstructed in Section 7,
  `Np = 2^(2p+2)·(p+1)·(2p+1)·pi^(2p+1) / (p·[(2p+1)!!]²)`, which was
  verified there against the explicit E1, E2 **and** E3 forms. This is where
  that reconstruction stops being unused and becomes load-bearing, so
  re-verify it against all six legacy values above before relying on it, and
  record the comparison.
- `EV_PER_CM` becomes `EV2CM` from the shared module. Note that the existing
  8.06554393734921 is already the CODATA value — it is right, and it is right
  by coincidence of maintenance, which is the problem.

Expect the six legacy constants to move by ~1e-5 to ~1e-6, in the same
direction and for the same reason as everything in Section 6; each shift must
be tabulated in the section against its legacy value, like every other Stage B
in this plan. `conv_out.py` and `calc_TP_M1E2.py` outputs move with them, and
`CODE/tests/test_conv_out.py` and `test_calc_tp.py` carry expected values that
will need updating — **check each changed expectation by hand before changing
it.** A test expectation edited to match new output is not evidence of
anything.

`CODE/perl_legacy/` is untouched. It is the reference for intent (CLAUDE.md
§8) and is not maintained; note in `radiate.py` that its constants are now
knowingly stale, so the next reader does not "restore" them.

Estimate: 3–4 h.

**What was done.**

- **`CODE/toolbox/fundconst.py`** — new, the recommended mechanism: a
  generated copy, not a parse. All 37 names of `FOR/fundconst.f`, values
  re-typed from the NIST listing, derivations written in the same
  association because float arithmetic is not associative and the tie is
  equality rather than a tolerance.
- **`CODE/tests/test_fundconst.py`** — new, six tests. It compiles
  `FOR/fundconst.f` into a driver that prints every parameter it declares in
  `1PD26.17` and asserts each equals its Python counterpart with `==`. **All
  37 agree bit for bit.** Two structural tests go with the values: a Fortran
  name with no Python counterpart fails, and a Python name that is not in the
  Fortran fails — so the copy cannot drift as the file grows and cannot
  become a second source of truth. The second is what will catch Section 13's
  four new entries. Skips without gfortran.
- **`CODE/toolbox/radiate.py`** — rewritten on top of it. E1, M1 and E2 are
  spelled exactly as `GAFROMSCOEF` spells them and equal the compiled
  functions **to the last bit**, asserted in both directions; that is the tie
  between the two languages the section asked for. The spelling matters: the
  Grant factor gives the E2 type factor as `5.52754142998615672D-17` against
  the Fortran's `…795D-17`, the same number algebraically and one bit apart.
- **M2, E3 and M3** take the general Grant factor. Only the *ratio*
  `Np(p)/Np(1)` is used, so the reconstruction's prefactor divides out and the
  absolute scale still comes from `GAFROMSCOEF`; for p=2 that ratio is
  `pi**2/5` exactly in floating point, and a test asserts it. Re-verified
  against all six legacy values before being relied on, as the section asked.
- **`CODE/tests/test_coefficients.py`** — seven tests appended, the toolbox
  layer.

| | | derived (CODATA 2022) | legacy | rel |
|---|---|---|---|---|
| gA E1 | λ³ | 2.026126884494784e+18 | 2.0261269e18 | −7.65e−09 |
| gA M1 | λ³ | 2.697350022055934e+13 | 269735e8 | +8.18e−09 |
| gA E2 | λ⁵ | 1.119950029645370e+18 | 1.11995e18 | +2.65e−08 |
| gA M2 | λ⁵ | 1.490971399808825e+13 | 14909714e6 | −1.28e−10 |
| gA E3 | λ⁷ | 3.144416603621083e+17 | 3.14441e17 | +2.10e−06 |
| gA M3 | λ⁷ | 4.186111077266150e+12 | 4.1861e12 | +2.65e−06 |
| gf E1 | λ¹ | 3.037556835275454e+02 | 303.75568885954 | −1.76e−08 |
| gf M1 | λ¹ | 4.043850392256858e−03 | 0.0040438504 | −1.91e−09 |
| gf E2 | λ³ | **1.679022125292271e+02** | **167.84224** | **+3.57e−04** |
| gf M2 | λ³ | 2.235255057986555e−03 | 0.002235255 | +2.59e−08 |
| gf E3 | λ⁵ | 4.714089833354381e+01 | 47.140897 | +2.83e−08 |
| gf M3 | λ⁵ | 6.275791715355832e−04 | 0.000627579 | +2.73e−07 |

The E1 gf shift is the one Section 11 said this section inherits.

**Defect found: the legacy `gf(E2)` is wrong, by +3.6e−04.** Not a stale
constant — an internal inconsistency in the legacy table. `gf = gA·λ²/(GACON·1e16)`
is exact and type-independent (`GFFROMGACOEF`, "all types"), so `1.11995e18`
and `167.84224` cannot both be right, and the Grant factor confirms the gA
one to 2.65e−08. `167.90221` is the value consistent with the rest of the
table. Found the same way as defect 0a.1 — by deriving what had been copied.
It reaches the E2 rows of `conv_out.py`'s absorption model, which move by
3.6e−04 against its E1 and M1 rows. Pinned by a test so the legacy number
cannot later be mistaken for a target.

**Defect found: `EV_PER_CM` was out by a factor of a thousand.** It held
8.06554393734921 — eV per *1000* cm⁻¹ — and both call sites in `conv_out.py`
divided by it a quantity the `*1000` had already converted to cm⁻¹. The Auger
`dE(eV)` column has been printing meV under an eV heading since it was
written. The Perl does the same (`conv_out.pl:733, 965, 968`) and its comment
names the constant "eV to cm-1", which *is* 8065.54: the value was the
mistake, the usage was right, in both languages. `WORK/AI/DIEL`'s own README
confirms it independently — its free electron is 18.58 Ry = 2038.92 kK =
**252.8 eV**, and the old expression gave 252 800. **Put to the author
2026-08-21, who chose to fix it**: both sites now take `fundconst.EV2CM`
directly and the column reads in eV. `radiate.EV_PER_CM` is *deleted* rather
than corrected in place, so a caller written against the old behaviour fails
loudly instead of quietly changing by 1000. A deliberate, documented
divergence from `CODE/perl_legacy/`.

**No test expectation needed changing, and that was checked rather than
assumed.** Every expectation in `test_conv_out.py` and `test_calc_tp.py` that
goes through `radiate` is a round trip, a ratio of two quantities converted
the same way, or a comparison of two rows of the same table — all insensitive
to the constant's value. Not one absolute S, gf or gA figure is written down
in either file. Nothing was edited, so nothing had to be justified.

Nothing in `FOR/` was touched and no executable was rebuilt.

Suite: **309 passed, 252 skipped** (289/252 before Section 11).

Record: `F:\COWAN\fundconst_work\py\README.txt`.


### Section 13 — The four remaining inline constants — **DONE 2026-08-21**

The author's standing rule, stated at the head of this plan: **every
fundamental constant comes from `fundconst.f`.** Section 10 found four that
do not, quantified them, and deferred them. This section closes that gap, and
with it the plan's scope statement becomes literally true rather than nearly
true.

| Site | Literal | What it is | CODATA 2022 | Legacy error |
|---|---|---|---|---|
| `rcg11kd.f` 1328 1329 (`SPRIN`), 8979 8980 9102 9103 (`CALCFC`), 15217 (`RCEINP`) | `0.50116D0`, `1.50116D0` | the electron g-factor as (g−1)/2 and 1+(g−1)/2, in the Landé g-factor | `0.501159652180`, `1.501159652180` | +6.9e-07, +2.3e-07 |
| `rcg11kd.f` 10454 (`MUPOLE`) | `1.00232D0` | g−1, in the magnetic-moment matrix element | `1.002319304361` | +6.9e-07 |
| `rcg11kd.f` 14917 (`RCOEFF`) | `2.1717D-8` | 2·√pi·alpha·c·a0², the Maxwellian rate coefficient per unit collision strength, cm³/s, T in rydbergs | `2.1716669877D-08` | +1.5e-05 |

Eight sites, four routines that need the include (`RCOEFF` already has it),
one program. Nothing outside `rcg11kd.f` is touched.

**New entries in `fundconst.f`.** One in Section A, four in Section B, in the
established form — CODATA value and uncertainty in the trailing comment, a
`C     LEGACY:` line above each Section B entry carrying the literal bit for
bit, derived quantities computed rather than transcribed.

- Section A: the electron g-factor, `-2.00231930436092` as published,
  entered as its magnitude with the sign convention stated in the comment —
  the code uses `|g|` throughout and a sign flipped here would be silent.
  u(r) = 1.3e-13, the most precisely known number in the file.
- Section B: `(|g|−1)/2`, `1 + (|g|−1)/2`, `|g|−1`, each derived from the
  Section A entry, and the Maxwellian coefficient
  `2·√pi·alpha·c·a0²` derived from `ALPHFS`, `CLIGHT` and `ABOHR`.
  Verified: 7.2973525643e-3 · 2.99792458e10 · (5.29177210544e-9)² ·
  2√pi = 2.1716669877e-08, against the 2.1717e-8 in the source.

The three g-factor names are three spellings of one constant, exactly as
`RINF` / `RINFKK` and `ALPH2H` / `ALPH2Q` already are; that is the file's
existing idiom and there is no reason to invent another. Candidate spellings
`GEFACT`, `GEM1H`, `GEP1H`, `GEM1` and `CMAXWL` — **not settled**, because
the collision check settles them. Re-run it over `FOR/*.F`, `FOR/*.f` and
`FOR/LSJJ/` before writing any of them down: this is the rule at the head of
`fundconst.f`, it has already forced two renames in this plan (`ALPHFS`,
`PICON` in Section 3; `TWOPIC` in Section 10), and a `PARAMETER` that shadows
a read-only local is a silent wrong answer, not a compile error.

**Includes.** `SPRIN`, `CALCFC`, `MUPOLE` and `RCEINP` need
`INCLUDE 'fundconst.f'` added after their `IMPLICIT` card. Check each for a
local of the chosen names — `CALCFC` and `RCEINP` are long routines with
large local name spaces — and check that none already includes the file
under a path that would double-declare it.

**Protocol: the standard two stages (§3), on five regression cases.**

| | | |
|---|---|---|
| `stageA` | includes added, eight sites substituted, parameters at LEGACY values | must be **byte-identical to the committed build** |
| `stageB` | values updated in `fundconst.f` only | the value change, attributable by construction |

Stage A here is worth more than usual, because adding an include to three
routines that never had one is the exact circumstance in which a shadowed
local goes unnoticed. Run all five — `SnVIII`, `Forbid`, `AI`, `AI/DIEL`,
`Pr3` — and require md5 identity on `TAPE2E`, `OUTGINE`, `TAPE19`,
`genealogy.txt` and on `OUTG11` less its timing lines.

**Stage B moves printed output, and this is the section's one visible
change.** Unlike Sections 4–6 the shift lands in a quantity the author reads
directly:

- The **Landé g-factors** printed by `SPRIN` (`WRITE (IW,20)`), written to
  `TAPE19` by `CALCFC` (`WRITE (19,91)`), and rebuilt in `RCEINP`, all rise
  by about 2.3e-07 relative — well below the printed field width, so the
  expectation is that **no printed g-factor moves at all**, and any that does
  is a knife-edge rounding in the last decimal. Verify, do not assume:
  `cmp_num.py --rtol 1e-12` over `OUTG11` and a direct diff of `TAPE19`.
- The **M1 matrix elements** through `MUPOLE` move by 6.9e-07 relative,
  carrying the M1 `gA` and `gf` with them. `Forbid` is the only deck that
  exercises M1 and is the case that decides this section.
- The **Maxwellian rate coefficients** through `RCOEFF` move by −1.5e-05,
  the largest of the four and the same order as the Section 6 shifts. This
  reaches whatever deck turns the collision-strength path on; if none of the
  five does, say so plainly in the section the way Section 10 said it of
  `PICSCN`, and check the value algebraically instead. Do not let an
  unexercised site be reported as a verified one.

Tabulate every shift against its legacy value, as every Stage B in this plan
has.

**Then delete the holding pen.** `PHYSICS_STILL_INLINE` in
`CODE/tests/test_literals.py` goes away entirely — the group, its docstring,
and its entry in the `INVENTORY` loop — and the four tokens move into
`CODATA`, which is the group for literals that live in `fundconst.f`. The
module docstring paragraph that explains the holding pen goes with it. After
this section the test asserts something stronger than it does today: that
`FOR/` contains **no** physical constant outside `fundconst.f`, with no
exception list. Add a sentence to the test's header saying exactly that, so
that a future reader who wants to add an exception has to argue with it.

**And say it in the source.** The `rcg11kd.f` header change list gets its
entry, and the header of `fundconst.f` gets one line stating the rule the
author set: a fundamental constant belongs in this file, and there is no
threshold of smallness below which it does not. Section 9 then documents the
completed state without qualification.

Estimate: 2–3 h, the bulk of it the five-case Stage A verification.

**What was done.** Work area `F:\COWAN\fundconst_work\gfac\`, full record
in its `README.txt`.

- `GEFACT = 2.00231930436092` in Section A of `fundconst.f`, and `GEM1H`,
  `GEP1H`, `GEM1`, `CMAXWL` in Section B, each derived from it (or, for
  `CMAXWL`, from `PICON`, `ALPHFS`, `CLIGHT` and `ABOHR`) rather than
  transcribed. The collision check came back clean on all five candidate
  names, so unlike Sections 3 and 10 no rename was forced and the spellings
  proposed above stand. `DSQRT` inside a `PARAMETER` is a gfortran
  extension; it was verified to compile and to fold before `CMAXWL` was
  written that way.
- Eight literals substituted in `rcg11kd.f`; `SPRIN`, `CALCFC`, `MUPOLE` and
  `RCEINP` given `INCLUDE 'fundconst.f'`. Each was checked against the whole
  name list of the file, not only the four new names. Zero shadowed locals.
- **RCN and RCN2 were not re-run, and did not need to be.** Both include
  `fundconst.f`, so the change could in principle reach them. It does not:
  `RCN36K.F` and `RCN2K.F` were compiled against the old and the new
  `fundconst.f` with the production flags and the **object files are
  bit-identical**. That is a stronger statement than a matching run, and it
  also removes the Norton problem of Section 10 — no RCN2K binary is built
  under the work area at all. The case directories are seeded with the
  Section 11 baseline's RCN and RCN2 output and only RCG is run.
- Stage A (`sweep/stageS11` vs `gfac/stageA`): **0 files differ** on all
  five cases. Stage B rebuilt from the final source after the header and
  comment edits (`stageB` vs `stageB2`): **0 files differ**, so the prose
  changed nothing.

**Result — the value change, field by field.** Seven files move in five
cases, and every one of them is accounted for.

| | legacy | CODATA 2022 | rel |
|---|---|---|---|
| `GEM1H` | `0.50116` | 0.50115965218046 | −6.94e-07 |
| `GEP1H` | `1.50116` | 1.50115965218046 | −2.32e-07 |
| `GEM1` | `1.00232` | 1.00231930436092 | −6.94e-07 |
| `CMAXWL` | `2.1717e-8` | 2.1716669877207e-08 | −1.52e-05 |

- **`OUTG11`, Sn VIII — one field in 470 kB.** The J=1.5 fourth level's
  Landé g-factor printed `2.100` and now prints `2.099`. That is the
  knife-edge rounding the section predicted, and it was verified rather than
  assumed: a throwaway build with `F20.14` in place of `F9.3` reads the
  number as 2.09950010910943 → 2.09949934602454, −3.6e-07, crossing 2.0995.
  It is the only printed g-factor that moves anywhere.
- **`OUTG11`, Forbid — two fields**, both `MUPOLE` sum-rule totals over
  squared matrix elements: `SUMS21` −2.30e-07 and `SUMS20` −1.35e-06, the
  latter bounded by 2 × 6.94e-07 because the M1 element enters squared. No
  `gf`, no `gA`, no M1 or E2 matrix element and no line count moved; the
  printed fields are too narrow to show 1.4e-06. Forbid is the deck that
  decides this section, and it decides it quietly.
- **`OUTG11`, Pr III / AI / AI DIEL** — unchanged.
- **`TAPE2E` — differs in all five.** The only Section 13 quantity in it is
  the g-matrix `RCEINP` packs onto unit 15. An instrumented build dumping
  that matrix at `E26.18` measures the whole of it: max |Δ| 1.62e-06
  (Sn VIII, Pr III), 1.16e-06 (Forbid), 6.96e-07 (AI), 1.39e-07 (AI DIEL),
  on elements of magnitude ≤ 3.34. `OUTGINE` is the same information
  printed and is **byte-identical everywhere**, so nothing RCE reads at
  printed precision moves.
- `ING11`, `TAPE2N`, `OUT36`, `OUT2`, `genealogy.txt` — byte-identical.

**Coverage — three sites the runs do not exercise, established positively
rather than inferred.** Each was found with an instrumented build that
produced no output at all, not by reading the absence of a diff.

- `CMAXWL`: `RCOEFF` is not reached by any deck in `WORK/`. Checked
  algebraically instead — 2√π·α·c·a0² = 2.1716669877207e-08 against the
  `2.1717e-8` the source carried. Unexercised, and said so, as Section 10
  said it of `PICSCN`.
- **`TAPE19` is never produced.** The section above asked for md5 identity
  on it; `CALCFC`'s `WRITE (19,91)` sits on a path none of the five decks
  reaches. Stage A proves the substitution there; the value change is
  unexercised.
- `CALCFC`'s `P5` g-factors, likewise: that block is under `IF (KCPL.LT.2)`
  and every deck here runs `KCPL=2`.

**Two things were repaired in passing.**

1. **`fundconst.f`'s site line numbers had gone stale.** Section 11 moved
   `rcg11kd.f` by about sixty lines without updating them, and this
   section's header entry moved it again. All 23 lists were recomputed from
   the current source and rewritten — `RCN36K.F` and `RCN2K.F` included,
   where they were out by 30 to 45 lines. Two facts surfaced while doing it
   and are now recorded at their entries: **`ALPHIN` has no consumer left in
   RCG at all** — its one site went with the inline `PSQ` derivation in
   Section 11, and the file still claimed 13188 — and `PICON` and `GFE1` are
   now used only inside `GAFROMSCOEF`.
2. `fundconst.f`'s header now carries the author's rule in as many words,
   and `CODE/tests/test_literals.py` enforces it with **no exception list**:
   `PHYSICS_STILL_INLINE` is gone — the group, its comment block, the
   docstring paragraph and its entry in the `INVENTORY` loop — `GEFACT`'s
   literal joins `CODATA`, and the four tokens join the un-substitution
   guard, now `test_the_substituted_literals_have_not_come_back` over
   `SECTION_10` plus `SECTION_13`.

**And the Python side stayed in step.** `CODE/toolbox/fundconst.py` gained
the same five names in the same order and the same association;
`test_fundconst.py` compiles the include file and compares bit for bit, and
`CMAXWL` agrees exactly, so `math.sqrt` and `DSQRT` round the same way here.
Full suite: **309 passed, 252 skipped.**

---

## 5. Effort

| Section | Estimate |
|---|---|
| 1 comparator tool | 2–3 h — **done** |
| 2 baselines | 1 h (mostly Pr III runtime) — **done** |
| 3 `fundconst.f` | 1 h — **done** |
| 4 RCN | 3–4 h (the single-precision hazards) — **done** |
| 5 RCN2 | 1 h — **done** |
| 6 RCG constants | 1.5 h — **done** |
| 7 RCG coefficient routines | 2–3 h — **done** |
| 8 RCE | 15 min — **done** |
| 10 sweep for missed sites | 2–3 h — **done** |
| 11 RCG onto the routines | 3–4 h — **done** |
| 12 Python side | 3–4 h — **done** |
| 13 the four remaining inline constants | 2–3 h — **done** |
| 9 docs (runs last) | 2 h (incl. §0c reminder x3) — **done** |
| final promotion of the four executables | 1 h — **done** |

**~3 working days**, roughly 25% editing and 75% verification. Sections 1–3
are prerequisites; 4 through 8 are strictly sequential, because each consumes
the previous stage's output. Sections 10, 11, 12 and 13 are independent of
one another and may be done in any order, but all four come before Section 9
and before the final promotion — 10, 11 and 13 all move `rcg11kd.f`, so
rebuilding the executables between them would only have to be undone.

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
6. **New, from Section 10: the autoionization rate is 0.05 % low today.**
   `2066` should be `2067.0687`. This is not a rounding refinement like the
   rest of this plan — it is sixty times the size of the typical Stage B
   shift, and it changes numbers the author publishes. The change itself is
   not in doubt; what needs the author's decision is the same question as
   open question 3, one level up: does the corrected autoionization rate
   become the benchmark from here on, and is anything already published
   against the old one? Nothing downstream is blocked on the answer, so
   Section 10 proceeds either way; only the closing note changes.
7. **New, from Section 11: preserve `SPECTR`'s 1e8 wavelength scaling in the
   E2 path, or remove it?** Keeping it means the rewired code still carries
   an introduced-and-cancelled factor of 1e16 and still differs from
   `GAFROMSCOEF(3)` in the last ulp; removing it makes the two identical and
   the arithmetic honest, at the cost of a wider edit to a routine nobody
   wants to disturb. Recommendation: **keep it** for Section 11 and note the
   ulp, because the section's value is in removing the *duplicate
   derivation*, and a wider edit puts the byte-identity acceptance test at
   risk for no numerical gain. Revisit only if a diff appears that the ulp
   explains.

8. ~~**New, from Section 10: four physical constants are still written inline
   — fold them in, or leave them?** Each is an order of magnitude or more
   below the `2066` error, and each moves printed output, including the
   printed Landé g-factors.~~ **Resolved (author, 2026-08-20): fold them in.**
   The rule is not a threshold — every fundamental constant comes from
   `fundconst.f`, whatever its error. Section 13 carries it out, and no part
   of this document may say otherwise. This also settles the same question
   for anything a future sweep turns up: the answer is already given.
