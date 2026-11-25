# Hydrogen Atom Simulation: An Early Undergraduate Exploration

This archive preserves an early personal exploration of quantum mechanics through study notes, small Python experiments, and visual media. It is a record of learning, not a finished scientific software package.

## A personal note

In the first semester of my Engineering Physics degree, we worked through the Schrödinger equation for a one-dimensional potential well. I found the subject fascinating and wanted to push further, so in January 2023 I worked through the three-dimensional version on my own using online sources. I deliberately stayed away from AI while doing this: I wanted the challenge of wrestling with the derivations myself.

Shortly afterwards, I wrote code to turn some of what I had learned into orbital-style visualizations. My laptop at the time could not conveniently render a video directly, so I used a simple workaround: generate many individual snapshots and combine them into a video. I also made a four-sided, symmetric layout so I could play it on my phone, place an inverted plastic pyramid over the screen, and view the orbitals as a small hologram.

Looking back, the physics and code have real shortcomings; this archive does not try to hide or silently repair them. What matters to me is that this was my first serious attempt to connect quantum theory, programming, and a physical visualization. It taught me how much care is needed when translating equations into code, how to work around limited hardware, and why an attractive visualization is not the same thing as a validated scientific result.

## What's included

| Folder or file | Contents |
| --- | --- |
| [`notes/handwritten/`](notes/handwritten/) | Your handwritten hydrogen-atom notes in PDF form. |
| [`notes/typed/`](notes/typed/) | Markdown transcription and organization of the derivations. |
| `notes/supplementary/` | An additional Word document retained locally and not included in this public release. |
| [`src/`](src/) | Copies of the original Python scripts and their dependency list. The code has not been corrected or modernized in this archive. |
| [`media/`](media/) | Selected original illustrations, the [hologram-pyramid template](media/hologram-pyramid-template.jpg), and the compact [four-panel hologram demo](media/hologram-pyramid-demo.mp4). |
| [`ERRATA.md`](ERRATA.md) | Known technical limitations and points that need checking. |

The original workspace contains many more generated frames and video exports. They are intentionally not duplicated here; only one representative four-panel video is included.

The handwritten PDF is included in the public archive. The supplementary Word document remains local and is ignored by Git.

## Reading the notes alongside the code

The typed notes cover the 3-D infinite potential box, angular momentum in spherical coordinates, and the hydrogen atom. The main orbital visualizer contains hand-written expressions for selected states, rather than a general implementation of the complete wavefunction derived in the notes. The separate finite-well script is a one-dimensional numerical experiment; it is not the 3-D box implementation.

## Running the archived scripts

[`requirements.txt`](requirements.txt) lists the Python packages used by the scripts. The orbital visualizer is preserved as an archival copy, not packaged as a turnkey application: it renders many frames by default and expects a state-specific output folder relative to its working directory. The generated frame folders are not included here. See [ERRATA](ERRATA.md) before interpreting the scripts' outputs.

## Before public redistribution

This is a personal learning archive, not verified physics software. No license is included.
