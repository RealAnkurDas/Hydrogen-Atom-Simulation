# Errata and known limitations

This file documents issues noticed while reviewing the archived material. It is not a corrected derivation or a scientific validation. The Python files are copied as they were found; no physics or behavior has been changed in this archive.

## Orbital visualizer

File: [`src/Hydrogen Electron Probability Cloud.py`](src/Hydrogen%20Electron%20Probability%20Cloud.py)

- **Bohr-radius scale differs from the notes.** The script uses `a_0 = 0.529e-11`, while the typed notes state `0.529e-10 m`. Those values differ by a factor of ten.
- **The n=3 exponential expressions have suspicious grouping.** The `300`, `310`, and `311` branches use `math.exp((-2*r)/3*a_0)`. As written, this multiplies by `a_0`; it does not divide by `3*a_0`. This is a likely implementation error, not a correction applied here.
- **The plotted quantity is not established as a normalized probability density.** Each state has a separate hand-written expression; some branches take square roots or repeated roots, and some values are negated. These transformations may make a visual effect, but no normalization or comparison against a trusted analytic implementation is performed.
- **Angular-coordinate handling is fragile.** Expressions such as `atan(y/(x + eps))` do not reliably preserve all quadrants, and epsilon offsets can perturb values. The code does not test coordinate edge cases.
- **The calculation is specialized and expensive.** The main block renders 360 views of state `321`. It builds a large Cartesian product before applying a radius filter and evaluates points in Python loops.
- **Output paths depend on the current working directory.** The script saves into a `Psi...Images` folder using a Windows-style relative path. The destination folder must exist, and running from a different directory or operating system may not work as expected.
- **Frames and videos are not scientific validation.** Their appearance does not establish that the encoded wavefunction, probability distribution, or physical scale is correct.

## Finite potential well

File: [`src/Finite Potential Well Solutions.py`](src/Finite%20Potential%20Well%20Solutions.py)

- The script searches for roots of transcendental equations, but calculates energy levels using a fixed `k = 1.466 / L` rather than the corresponding root. The printed energies therefore do not follow from the roots it found.
- `fsolve` is started from a small set of guesses near equations with tangent/cotangent singularities. The script does not robustly constrain roots to valid intervals or reliably deduplicate numerically close solutions.
- The final additional energy expression is not documented or checked against the preceding derivation.

Treat its output as an experiment, not as dependable energy-level data.

## Notes and scope

- The Markdown notes preserve a learning derivation and source notation. They have not been independently proofread line-by-line for algebra, units, boundary conditions, or convention-dependent normalization factors.
- There is no general solver for the 3-D infinite box, no reusable spherical-harmonic implementation, and no automated test suite in this archive.
- `src/Hydrogen.py` is only a small plotting stub; it is not another hydrogen simulation.
- The archive includes one representative four-panel hologram video, not all original frame sequences and exports. The full output collection was much larger and some exports appear to be alternate versions.

## What this exploration taught me

- Deriving an equation and implementing it are separate tasks; units, parentheses, coordinate conventions, normalization, and numerical methods all need their own checks.
- A compelling image can still be mathematically or physically misleading.
- Generating frames separately and assembling them later was a practical response to limited hardware, and thinking about the display geometry led to the four-sided pyramid layout.
- Keeping the handwritten notes, typed derivation, code, and media together makes it easier to trace how the idea developed—and to see clearly where understanding was still incomplete.
