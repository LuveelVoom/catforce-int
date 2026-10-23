# CatForce INT

A modified **CatForce**, the catalyst/conduit search program for Conway's Game of Life by Michael Simkin,
extended to isotropic non-totalistic (INT/Hensel) rules, other Life-like rules and Golly `MAP` rules.
Original CatForce is Michael Simkin's work; this repo holds LuveelVoom's extended version.

Main additions: per-rule compiled SIMD kernels, Apple-silicon support, bait objects, search/exclude areas,
`also-required` cells, unique-outcome output, lightcone pruning, contact filter, pause/resume, and an
in-browser `.in` file generator. Every input keyword is documented in [CatForce/README.md](CatForce/README.md);
development notes and measurements are in [HANDOFF.md](HANDOFF.md).

| Path | What it is |
|---|---|
| `CatForce/` | Source. `cd CatForce && make`, then `./CatForce file.in` |
| `CatForce/README.md` | Input-file reference (every keyword) |
| `generator/catforce_in_generator.html` | `.in` file generator; open in a browser |
| `generator/src/` | Generator sources; `python3 generator/src/build.py` rebuilds the page (both copies) |
| `tests/` | Regression suite (checks a build against recorded results) |
| `diag/` | Python helpers to decode and simulate RLE results |
| `user_files/` | Example search files |

## Quick check

```
cd CatForce && make
./CatForce ../user_files/highlife2.txt --max-time 12     # ~1 s, 7 results
bash ../tests/suite.sh ./CatForce check                   # ends with "all match tests/base.txt"
```
