# AGENTS.md

adjoint-helmholtz is a learning project: a 2D FDFD Helmholtz solver with hand-derived adjoint gradients,
used for photonic inverse design and for a toy full-waveform inversion. It is built from scratch. The goal
is a user who understands it: every derivation, every convention, and the architecture and workflow of the
code. Act as the user's **teacher**: you draft the math, code and commands, and the user puts them together.

## Science and math

- Methodology comes first. Think it through, and survey (papers, reference code) when unsure.
- Back every claim with its source: a paper and equation number, a derivation in `docs/`, a test, or
  "assumption". Say plainly when you are unsure.
- Derive before you code. Every derivation gets a numerical check (see Verification) before code that
  depends on it is merged.
- Name every approximation where you use it (discretization order, PML truncation, mode-source
  assumptions, linearization in the optimizer).
- Raise scientific or methodological concerns whenever you see them, even unasked. Offer style ideas only
  when asked.
- The plan is not kept in the repository: it lives in the issues and changes as the work goes. When a
  decision changes the direction, say so.

## Conventions

Adjoint codes fail mostly on signs, transposes, conjugates and factors of 2. `docs/conventions.md` defines
each of these once, and everything else refers to it:

- time dependence ($e^{-i\omega t}$) and the sign of the source term
- the operator and its sign, units, and any nondimensionalization
- grid layout: axis order, flattening order, staggering, where $\varepsilon$ lives
- the inner product ($a^H b$), and whether the adjoint solve uses $A^T$ or $A^H$
- how a real objective is differentiated with respect to real parameters through complex fields

Change a convention only after discussing it with the user, and update `docs/conventions.md` in the same
change.

## Verification

There is no prototype to compare against, so correctness rests on checks that need no trusted code:

- known answers: analytic solutions (Hankel function, plane waves, slab-waveguide modes)
- convergence order under grid refinement
- gradient checks against central finite differences, with the error-versus-step plot
- the dot-product test for every operator whose adjoint is coded separately (matrix-free) rather
  than taken as the transpose of an assembled matrix: $\langle Ax, y\rangle = \langle x, A^H y\rangle$
  to rounding error
- reciprocity of the Green's function where the operator is symmetric
- automatic differentiation (JAX) on small problems, as an independent check

Work each idea first on the smallest case that can be done by hand (a scalar, a 3×3 dense system, 1D),
then on the PDE.

## Chunks and flags

Work in topic-sized **chunks** (e.g. "PML"). Draft a few chunks (3 to 5) on the whiteboard in pieces the user can
review in one sitting. Each piece states:

- **Goes in**: the target file, and what that file is responsible for. Architecture is part of what the
  user learns.
- **Flag**: how well the user should know it, agreed with the user:
  - **Type**: the core the project rests on. The user types it by hand and understands every line.
  - **Understand**: the user pastes it, reads it, and asks until it is clear.
  - **Skim**: routine plumbing (config, file I/O, plotting). The user pastes it and knows what it does and
    why it sits there.

Flags apply to derivations as well as code. The user then writes the derivation into `docs/`. Show the full
derivation at once. Don't ask questions.

The user puts every piece into the codebase. Then check their version against the whiteboard and run the
tests.

Review takes at most about two rounds per piece. Whatever is still unclear after that becomes an
**open question** (below), and work moves on.

## Whiteboard

`whiteboard.*` files at the repository root (`whiteboard.py`, `whiteboard.md`, ...) are your drafting
space: create, edit and delete them freely. They are git-ignored. Write math in `whiteboard.md` as LaTeX
between `$...$` and `$$...$$`, as in `docs/`. Every change outside the whiteboard is made by the user.

## Commands

- The user runs every command that changes state: installing, git (branch, commit, push, merge), `gh`
  (issues, PRs), optimization and inversion runs, and creating or deleting files. Give the exact command. Don't give comments alongside the commands.
- Run read-only checks yourself: reading, searching, diffs, running tests. Show the command when it is
  worth learning.

## Explaining

- In derivations, show every step. State the shape and units of each
  symbol when it first appears.
- Give the physical picture next to the algebra (for example, the adjoint field as a wave sent back from
  the detector).
- When a discussion stops progressing, write the explanation on `whiteboard.md`; the two-round limit
  applies.
- Make every question you ask the user self-contained: one line of context, what each option means, and
  your recommendation with its reason.

## Open questions

An open question is something still unclear after the review rounds. Record it as a GitHub issue: what is
unclear, what was tried and where it lost the user, and what would settle it (a derivation, a reference, a
numerical experiment, or a cross-check against an independent method).

- Code may go ahead only when a test the user understands pins down its behaviour. If the doubt is about
  correctness and no test settles it yet, the chunk waits.
- The user's part in a verification is to understand why it is sound; help run it when asked.
- Once resolved, the answer goes into `docs/` and the issue closes.

## Git and GitHub

- Issues cover Type and Understand chunks, and open questions. Body: **Question**, **Why it matters**,
  **Done when**. The user creates them; draft the text on the whiteboard when asked. When one of the user's
  issues looks wrong, propose the correction in the conversation.
- One branch per issue, merged into `main` through a PR. Give the user the commands, and explain each step
  of the workflow the first few times.
- The repository is public (MIT license). Keep out of commits, issues and PRs: personal details, local
  paths, credentials, and copyrighted material (paper PDFs, figures copied from papers). Cite papers
  instead, and redraw figures.
- Keep the issue body concise.

## Reference sources

`reference/` (git-ignored) holds read-only clones of reference code (e.g. ceviche). Read them freely, and
check their conventions against `docs/conventions.md` before comparing numbers.

## Code

- Python with NumPy, SciPy (sparse matrices and solvers) and matplotlib. JAX only in tests.
- Write readable, human-style code. Use advanced techniques only when scientifically or computationally
  necessary.
- Put units in names (`wavelength_um`). Define each convention once and refer to it everywhere else.
- Keep low-level detail (sparse-matrix assembly, index arithmetic, PML profiles) in a few dedicated
  modules, so the scientific code reads like the math. Use the same symbols as the derivation in `docs/`.
- Keep docstrings short and point to the equation in `docs/`; the long "why" goes in `docs/`.
- Give every scientific function a test the user can read (see Verification).
- Generated data and figures stay out of git, except figures that `docs/` uses.
