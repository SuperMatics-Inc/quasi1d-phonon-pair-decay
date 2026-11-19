# Pair density matrices and axial helicity for quantum bond phonons on cylinders

Davis Rens, SuperMatics Inc., ORCID 0000-0001-5839-5197.

Article: https://doi.org/10.5281/zenodo.24369897  
Supplement: https://doi.org/10.5281/zenodo.24369899  
Source: https://github.com/SuperMatics-Inc/quasi1d-phonon-pair-decay

The note gives bounds for positive-temperature cylinders and chains with quantum bond phonons. The Holstein result is of known type. The work does not establish a square-lattice BKT bound. `RIGOR.md` contains the extended proof and the source is in `note/`.

## Exact checks

With Python 3:

```sh
python -B replay.py --all-mutants
python -B replay_pressure.py --all-mutants
```

The two aggregates reject 10 cylinder and 9 pressure mutants. These exact arithmetic checks cover selected finite coefficients and logic, not the analytic proof.

## Lean

With elan and the pinned Lean toolchain:

```sh
cd lean
lake exe cache get
lake build
lake env lean Core.lean
lake env lean Pressure.lean
```

These files formalize finite lemmas only, not the functional-analytic or thermodynamic arguments. Mathlib is pinned in `lean/lake-manifest.json`.

## Manuscript

With `pdflatex` and BibTeX installed:

```sh
cd note
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

AI research agents generated proofs, checks, code and manuscript drafts under D. Rens' direction and review. The software is released under Apache-2.0; the manuscript is CC BY 4.0.
