# Files in the distribution

[Back to the README](../README.md)

1) Root directory        - The front page and the terms of use.

	Contents:

        ```
        README.md          - what the package is, quick start, and the
                             map of this documentation;
        LICENSE.md         - the terms the package is distributed under;
        NOTICE.md          - attributions and the provenance of each part.
        ```

2) DOCS directory        - The documentation.

	Contents:

        ```
        README.md          - index of the documents below;
        INSTALL.md         - requirements and installation, all systems;
        BUILDING.md        - recompiling the Fortran, and other compilers
                             and operating systems;
        MANUAL_ADDITIONS.md - what differs from Cowan's manuals: running
                             the programs, file formats, new options;
        UTILITIES.md       - the Python utility codes;
        CHANGES.md         - release notes, program versions, the
                             fundamental constants, the bugs fixed;
        FILES.md           - this list.
        Cowan/             - R. D. Cowan's own documentation, copied from
                             his directory at LANL in December 1994 and
                             preserved unaltered:
            RCN_DOC.txt      - manual for the RCN, RCN2 and HF8 programs;
            RCG_DOC.txt      - manual for the RCG program;
            RCE_DOC.txt      - manual for the RCE program;
            readme.cowan.htm - R. D. Cowan's write-up for his own version
                               of the package, December 1999.
        ```
 
3) CODE directory   - the executable files, the data files they need,
                       and the utility codes.

	Contents:

	a) Main executables:

        ```
        RCN36K.exe   - Single-configuration radial wavefunctions
        RCN2K.exe    - Single-configuration and configuration-interaction
                       radial integrals and transition matrix elements
        rcg11kd.exe  - Angular matrix elements, energy levels, spectra
                       and cross sections
        RCE20KD.exe  - The least-squares fitting program
        ```

        The `.exe` files are the Windows builds. `CODE/unix/build.sh`
        puts the Linux and macOS builds beside them under the same
        names without the `.exe`.
	
	b) Main wrappers for running the main executables:
 
        ```
        rcn.bat       - batch file for running the RCN program;
        rcn2.bat      - batch file for running the RCN2 program;
        rcg.bat       - batch file for running the RCG program;
        rce.bat       - batch file for running the RCE program;
        make_cfp.bat  - batch file for creating the binary data decks.
        mkverres.bat  - called by the two build scripts; compiles the
                        Windows version resource from version.txt.
        mark_time.bat - called by the four above; times the run.
        update11.bat  - the post-fit cycle, described in UTILITIES.md.
        ```

        The Linux and macOS counterparts are in `CODE/unix/`, one per
        `.bat` file and named the same without the extension, together
        with `build.sh` (compile), `install.sh` (install and set the
        permissions) and `_common.sh` (the banner and the timing, which
        is what `mark_time.bat` does on Windows).
 
	c) Data files common to all calculations:

        ```
        SENIOR        - Text file used by RCG.
        version.txt   - The package version and the identity strings the
                        Windows build compiles into the executables as a
                        VERSIONINFO resource. Read by mkverres.bat and
                        by nothing that runs. install.sh installs it on
                        Linux and macOS, where it is the only record of
                        the version.
        ```

        RCG reads four more files from this directory, and none of them
        is distributed: the binary CFP decks `FOR072`, `FOR073` and
        `FOR074`, and `LSJJ.DAT`. They are built on your own machine by
        `make_cfp` -- step 6 of the Windows installation procedure, and
        on Linux and macOS a thing `install.sh` does for you.

        `LSJJ.DAT` is the binary deck of Cowan's LS coefficients of
        fractional parentage together with the jj ones, from which RCG
        builds the true LS -> jj transformation. It is read only by a run
        that asks for true jj coupling (KCPLD(8) = 8 in column 38 of the
        ING11 main card), so a missing file is not an error until then --
        which is what makes it safe for its generator, `CODE/lsjj/`, to
        need `numpy`, the one optional dependency in the package.

        Once built, the three CFP decks are portable as bytes between
        gfortran targets, so a deck built on one of the three supported
        systems is read unchanged on the others; they have to be built
        again only for an `RCG` compiled with a compiler whose
        unformatted record layout differs from gfortran's. `LSJJ.DAT` is
        not portable in that way -- a different LAPACK moves a few dozen
        of its coefficients in the last bit, invisibly to everything RCG
        prints -- which is a further reason it is built where it is read.

	d) Compilation scripts and the gfortran runtime libraries:

        ```
        gfortran_O3.bat     - Batch file for compiling a source file with
                              the GNU gfortran compiler. Usage:
                              gfortran_O3 <source_file_name>
                              (The name must be given without the extension,
                              which must be .F)
        gfortran_O3_OMP.bat - Same, with OpenMP enabled. RCG and RCE are
                              threaded; the result is bit-identical to the
                              serial build at any thread count, because no
                              parallel region reassociates a sum.
        gfortran_dbg.bat    - Same as gfortran_O3, but with options for
                              debugging with gdb (-O0 -ggdb3, no OpenMP).
        mkverres.bat        - Called by the first two, not by the third.
                              It writes a VERSIONINFO resource from
                              `version.txt` and compiles it with `windres`
                              (which ships with the same MinGW as
                              gfortran); the build links the result. Any
                              failure there is harmless -- it prints one
                              line, produces nothing, and the build goes
                              on without the resource. There is no Unix
                              counterpart: ELF has no VERSIONINFO.
        libgfortran-5.dll   | The gfortran runtime, as distributed by
        libgcc_s_seh-1.dll  | WinLibs 13.2.0 r8. All six are needed:
        libquadmath-0.dll   | RCG and RCE are OpenMP builds and import
        libwinpthread-1.dll | libgomp, which in turn imports libdl.
        libgomp-1.dll       |
        libdl.dll           |
        ```
		(These are the libraries the distributed executables were linked
		against. If you rebuild with a different compiler, replace all six
		with that compiler's own, and do not mix the two sets. They are
		Windows files, and are not installed on Linux or macOS: there
		libgfortran and libgomp come from the system.)

		`CODE/unix/build.sh` is the Linux and macOS counterpart of the
		three batch files above, with `--serial` for a build without
		OpenMP and `--debug` for one with `-O0 -ggdb3`. Unlike the batch
		files it takes any number of source stems, and none means all
		four.

	e) Auxiliary files:

        `ING11.CFP`     - the complete CFP input deck for `RCG`. This is
                          what `make_cfp` feeds to `RCG` to build the
                          binary CFP decks `FOR072`, `FOR073` and
                          `FOR074` at installation, and to `mkdeck.py` to
                          build `LSJJ.DAT`; it is needed then and rarely
                          afterwards.

        `version.txt`   - the package version, the one place it is
                          written down, and the identity strings the
                          Windows build compiles into the four
                          executables. Nothing at run time reads it.

	f) Utility files: described in [UTILITIES.md](UTILITIES.md).

        g) unix/ - the Linux and macOS wrappers:

                       One file per command, carrying no extension --
                       they are typed as `rcg`, `add_eav`, `conf_list`
                       and so on, exactly like the `.bat` files without
                       theirs. Each of the four main wrappers and each
                       of the tool wrappers is the counterpart of the
                       `.bat` file of the same name, and does the same
                       thing. Three files are not counterparts of
                       anything on Windows:

                       build.sh    - compiles all four programs. Also
                                     build_serial.sh (no OpenMP) and
                                     build_debug.sh (-O0 -ggdb3); the
                                     switches --serial and --debug on
                                     build.sh do the same.
                       install.sh  - copies the package to a prefix and
                                     sets the permissions that make the
                                     read-only code directory a
                                     guarantee rather than a
                                     convention. --user, --system,
                                     --prefix DIR, --dry-run.
                       _common.sh  - sourced by the four main wrappers;
                                     the banner and the timing, which
                                     is mark_time.bat's job on Windows.

                       Every one of them carries its reasoning in a
                       comment header, which is the detailed
                       documentation of what it does and why.

4) FOR directory    - Fortran source files of the Cowan programs.
                       Fixed-form Fortran 77, built with -std=legacy.

        RCN36K.F
        RCN2K.F
        rcg11kd.f
        rcgpar.f      - Include file for rcg11kd.f (array dimensions)
        fundconst.f   - Include file for RCN36K.F, RCN2K.F and
                        rcg11kd.f: every fundamental constant the
                        package uses, as named PARAMETERs, at CODATA
                        2022. See "The fundamental constants" in
                        DOCS/CHANGES.md.
        RCE20KD.F
        RCE20KD_PAR.F - Include file for RCE20KD.F (array dimensions)
        LSJJCOM.F     | Include files for the true jj coupling
        LSJJID.F      | added to RCG and RCE in 2026
        LSJJPAR.F     |
        RCEJJ2E.F     - Include file for RCE20KD.F, true jj coupling
        RCEREL.F      - Include file for RCE20KD.F: the state of the
                        relative-energy groups (see DOCS/MANUAL_ADDITIONS.md).
        LSJJ/         - The record of the LS -> jj derivation: the
                        scripts that validated the coefficients against
                        Cowan's own, their logs, and the JJCFP deck. The
                        generator itself lives in CODE/lsjj/, where
                        make_cfp can reach it in an installed tree. Its
                        README.txt is the long-form account of that work,
                        and of everything else done to these sources
                        since: the utility-code port, the constants work
                        described in DOCS/CHANGES.md,
                        and the Windows version resource.


5) WORK directory   - Sample directory for Cowan calculation projects,
                       and the examples that come with the package. The
                       two files in WORK itself are the installation
                       test; each sub-directory is a worked example of
                       one thing the codes can do.

        IN36          - sample input file for RCN (for two configurations
                       of Sn VIII);
        IN2           - sample input file for RCN2;

        AI/           - B IV example of an autoionization calculation:
                       IN36 and IN2 only, run the four programs in an
                       empty directory holding a copy of them.

        AI/DIEL/      - the same physics set up by RCN's DIEL option,
                       which writes the continuum configuration cards
                       for you instead of your punching them by hand.
                       Inputs (indiel, in36, in2, gcard) and reference
                       outputs to compare your run against. See its
                       README.txt.

        Forbid/       - Sn VIII 4d7, forbidden lines: an RCG run with
                       the M1 and E2 spectra switched on (IMAG and
                       IQUAD, columns 49-50 of the IN2 main card), and
                       the example for calc_TP. IN36 and IN2 only.
                       Its README.txt explains why DMIN, columns 61-65,
                       must be 0.00 for forbidden transitions -- with
                       the usual electric-dipole cut of 5.00 or more,
                       RCG prints the section headings and no lines
                       under them.

        Mo15/         - Ni-like Mo XV: the relative-energy group example.
                       Twenty-two known levels, of which the five
                       3d9.4f and 3d9.5f ones are known to 400-600 cm-1
                       while the 3d9.4s and 3d9.4p ones are known to
                       under 0.1 cm-1, and ZETA(3d) is linked across
                       every odd configuration -- so the loosely tied
                       levels pull a parameter the well-measured ones
                       depend on. IN36, IN2, the experimental level
                       list, and RCEINP.example, the finished deck with
                       its letters and group cards. Its README.txt
                       gives the run step by step and compares the fit
                       with and without the weighting; it is also the
                       guidance on choosing SIGMA_INT.

        Pr3/          - Pr III, a real least-squares fit, and the worked
                       example of everything downstream of one: the
                       RCEOUT/RCEINP input format, the update11 cycle,
                       the merge of a coupling computed in a second
                       folder (Pr3/JK, the same fit asked for in
                       LS + LSJK instead of LS + JJ), and the three
                       programs of section f.2. It is also the conv_out
                       and calc_TP example: Pr3/JK carries a worked
                       calc_TP.par, and the README gives the two deck
                       fields an uncertainty run wants changed. Inputs
                       only -- the run writes about 412 MB, and a
                       calc_TP run of 25 trials about 600 MB more. It
                       also carries the dlv.dat and numset.dat an IDEN
                       run needs. Its README.txt gives the command
                       sequence step by step, with the figures each
                       step should report.

        W70/          - W LXX, boron-like tungsten: the true jj example
                       (KCPLD(8) = 8, column 38 of the ING11 main
                       card), and with it the true-jj_keys.txt key to
                       the 8-character basis labels and the LEVELS*J
                       files. It is also the example of a run with
                       nothing observed, so there is no fit, every
                       parameter is fixed, and make_levels_list_pub
                       has to be told --nlev. IN2 and IN36 only; the
                       run takes seconds. See its README.txt.

        U90/          - U XC, lithium-like uranium: a small, fast
                       least-squares fit over the K-shell-excited
                       configurations, and the reorder_ing11 example.
                       U90/reordered is the same fit with the ING11
                       subshell slots permuted, and shows what a
                       different coupling order does to the level
                       designations and to the percentages -- the same
                       eigenstates and the same energies, expanded in
                       a different basis. Inputs only. See its
                       README.txt.


                                  August 26, 2026

                                  Alexander Kramida,
                                  e-mail akramida @ umd.edu

