# License

This package is not covered by a single licence. It is a composite of four
strata with different origins and different terms, and the terms of each
stratum still apply to it:

| Component | Terms | Section |
|---|---|---|
| The original Cowan programs (RCN, RCN2, RCG, RCE) | 1980 Los Alamos / U.S. Government notice — free distribution, **no sale** | [1](#1-the-original-cowan-programs) |
| Portions developed at NIST (through the 2021 release) | NIST software notice; not subject to copyright within the United States | [2](#2-portions-developed-at-nist) |
| Work by A. Kramida from 2022 onward | MIT License | [3](#3-new-work-2022-2026) |
| Third-party components | see each item | [4](#4-third-party-components) |

**The practical effect** is stated in [section 5](#5-what-this-means-in-practice).
In short: you may use, copy, modify and redistribute this package and works
derived from it, anywhere in the world, at no charge — but you may **not sell**
it, because the programs at its core have carried that restriction since 1980
and this package cannot relicense them.

Because of that restriction the repository as a whole is **not** MIT-licensed
and should not be labelled as such. Only the portions described in section 3
are offered under the MIT License.

---

## 1. The original Cowan programs

`FOR/RCN36K.F`, `FOR/RCN2K.F`, `FOR/rcg11kd.f` and `FOR/RCE20KD.F` descend from
the atomic-structure programs written by Robert D. Cowan at Los Alamos. Each
retains its original notice in its header. `FOR/rcg11kd.f` carries the fullest
form:

> PROGRAM COPYRIGHT 1980, LOS ALAMOS NATIONAL LABORATORY AND U.S. GOVERNMENT.
> THIS PROGRAM MAY BE FREELY DISTRIBUTED GRATIS ANY PLACE IN THE WORLD, BUT MAY
> NOT BE SOLD.

The other three carry the same restriction in slightly different words. Those
notices are to be kept intact in any copy or derivative, and the prohibition on
sale travels with the code.

## 2. Portions developed at NIST

Substantial parts of this software were developed by Alexander Kramida while
employed at the National Institute of Standards and Technology, and were
distributed by NIST through the 2021 release
(<https://doi.org/10.18434/T4/1502500>). That work includes, among other things,
the packed binary file format (`PK`/`UNPK`), eigenvector recognition in RCE, the
Perl utilities now retired to `CODE/perl_legacy/`, and a large number of
corrections to the original programs. Software developed by NIST employees is
not subject to copyright protection within the United States. Those portions
were distributed, and continue to be distributed, under the following notice:

> NIST-developed software is provided by NIST as a public service. You may use,
> copy, and distribute copies of the software in any medium, provided that you
> keep intact this entire notice. You may improve, modify, and create derivative
> works of the software or any portion of the software, and you may copy and
> distribute such modifications or works. Modified works should carry a notice
> stating that you changed the software and should note the date and nature of
> any such change. Please explicitly acknowledge the National Institute of
> Standards and Technology as the source of the software.
>
> NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY
> OF ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW,
> INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST
> NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE
> UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES
> NOT WARRANT OR MAKE ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR
> THE RESULTS THEREOF, INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY,
> RELIABILITY, OR USEFULNESS OF THE SOFTWARE.
>
> You are solely responsible for determining the appropriateness of using and
> distributing the software and you assume all risks associated with its use,
> including but not limited to the risks and costs of program errors, compliance
> with applicable laws, damage to or loss of data, programs or equipment, and
> the unavailability or interruption of operation. This software is not intended
> to be used in any situation where a failure could cause risk of injury or
> damage to property. The software developed by NIST employees is not subject to
> copyright protection within the United States.

## 3. New work, 2022-2026

Work carried out after the 2021 NIST release — the replacement of `EQUIVALENCE`
and static storage by dynamic memory management, true *jj* coupling in RCG and
RCE, the speed and memory work, OpenMP parallelization, the per-configuration
term cutoff, the Python utilities in `CODE/` and `CODE/toolbox/` with their test
suite, the port to Linux and macOS in `CODE/unix/`, and the fundamental-constants
work — is the work of Alexander Kramida at the Department of Astronomy,
University of Maryland, and is offered under the MIT License:

> Copyright (c) 2022-2026 Alexander Kramida
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

The MIT grant above, including its permission to sell, applies to these portions
taken on their own. It does not and cannot lift the restrictions of sections 1
and 2 from the files those sections govern, and the four main Fortran programs
contain all three strata inseparably.

## 4. Third-party components

**GRASP memory-management routines.** `ALLOC`, `RALLOC` and `DALLOC`, present in
both `FOR/rcg11kd.f` and `FOR/RCE20KD.F`, originated in development versions of
the GRASP atomic-structure software and retain their original authorship
comments. See `NOTICE.md`.

**R. D. Cowan's documentation.** Everything in `DOCS/` is a verbatim copy of
Cowan's original documentation, redistributed unaltered for historical purposes.
It is not covered by section 3 and is not to be modified.

**GCC and MinGW-w64 runtime libraries.** The Windows distribution ships
`libgfortran-5.dll`, `libgcc_s_seh-1.dll`, `libgomp-1.dll`, `libquadmath-0.dll`,
`libdl.dll` and `libwinpthread-1.dll` beside the executables, as built and
distributed by WinLibs. The first five are components of the GNU Compiler
Collection, licensed under the GNU General Public License version 3 or later
**with the GCC Runtime Library Exception**, which is what permits their
redistribution alongside software under these terms; `libwinpthread-1.dll` is
part of MinGW-w64 and is distributed under its own permissive terms. Neither
their licences nor their copyright holders are altered by this package, and
their source is available from the GCC and MinGW-w64 projects. See
<https://gcc.gnu.org/onlinedocs/libstdc++/manual/license.html> and
<https://www.mingw-w64.org/>.

## 5. What this means in practice

You may, at no charge and anywhere in the world:

- run this software for any purpose, including commercial research;
- copy it and redistribute it;
- modify it and redistribute your modifications, provided you keep the notices
  in sections 1 and 2 intact and state that you changed the software, with the
  date and nature of the change;
- publish results obtained with it, with the acknowledgment requested in
  `NOTICE.md`.

You may **not** sell this software, or sell a derivative work that contains any
part of the programs described in section 1.

See `NOTICE.md` for the provenance of the package and for how to acknowledge it.
