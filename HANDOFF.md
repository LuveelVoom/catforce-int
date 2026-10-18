# CatForce INT — handoff for a new chat

Upload this whole zip to the new chat and say: "Read HANDOFF.md first, then continue."
Everything needed to keep working is in this archive. Nothing depends on the old chat.

## 1. What this project is

A modified **CatForce** (catalyst/conduit search program for Conway's Life, by Michael Simkin)
extended to **isotropic non-totalistic (INT/Hensel) rules**, Life-like rules and Golly `MAP` rules.
The user runs searches mainly in rules such as `B3/S23-q4i`, `B3/S23-a5` and HighLife `B36/S23`.

Folder map of this zip:

| Path | What it is |
|---|---|
| `CatForce/` | Current source tree. Build with `make`, run `./CatForce file.in` |
| `CatForce/README.md` | Full input-file reference (every keyword) |
| `generator/catforce_in_generator.html` | Latest .in file generator (single HTML page, open in a browser) |
| `generator/src/` | Generator source pieces and the build script (see §5) |
| `user_files/` | The user's own .in files, results and attachments from the chat |
| `tests/` | Regression suite (`suite.sh`), its inputs and baseline results |
| `diag/` | Small Python helpers used to decode/simulate RLE results (`lib.py`) |
| `reference/LC1.com.html` | Forum post (Lightcone author) used as inspiration for speedups |

## 2. User preferences (keep following these)

- Be concise; never repeat points; short answers even on long-thinking models.
- No deep-sounding unfounded generalisations.
- Don't agree for the sake of agreeing; evaluate corrections on evidence.
- Treat each chat fresh (no context from other chats) — this handoff is the deliberate exception.
- State confidence (high / moderate / low) and the type of source.
- Recently: **"use as little usage/tokens as possible"**, minimal code changes, stop and report
  when something gets expensive, minimal testing unless asked.
- Don't update the generator unless asked (or a major input-format rework requires it).

## 3. CatForce: features added in this chat (all in `CatForce/`)

Input-file keywords (all documented in README.md):

- `rule <rulestring>` — Hensel/INT, B/S numeric, or Golly `MAP…` (B0 rejected).
- `bait rle x y [survive n] [absent n]` — object that must be disturbed, reappear and survive (≤ 8).
- `search-area x y w h` (repeatable, union) and `exclude-area x y w h`.
- `also-required rle x y` — cells that must stay alive every generation.
  **Fixed late in the chat:** these cells are now automatically added to the starting pattern
  (before, a `pat` without them found nothing).
- `keep-duplicates`, `keep-interacting`, `output-gap`, `threads` (0 = all cores), `max-time`,
  `max-results`, `output-list`, `include`, `part K/N`, `verbose`.
- `unique-output file` — third output: one square per distinct outcome
  (default `<output>_unique.rle`). Each square is cropped to that outcome's starting pattern.
- `lightcone no` — disable lightcone pruning (on by default; measured roughly neutral to
  slightly slower on 2-catalyst tests, i.e. 7–30% overhead; meant for 3+ catalysts).
- `lightcone-tight yes` — UNSAFE faster variant: aim lightcone only at mismatched cells at failure time.
- `contact-filter yes` — automatic contact filter: skip a placement if the new catalyst, with no
  further catalysts, never comes back within its max-active. **Best measured speedup:**
  190 s → 10.5 s with identical 8 results (3q4i-3 cut to 2 catalysts, max-gen 120).
  Can lose results where a later catalyst rescues an earlier one.
- `resume-from n` — pause/resume. On Ctrl-C or max-time the program prints `resume-from n` and writes
  `<input>.resume`. A resumed run writes outputs with suffix `_from<n>`.
- `spaceship-cutoff n|off` (default 10) — no catalyst placements triggered by escaped spaceships
  more than n cells from the rest of the reaction; positions they pass are masked.
- `torus 64|128` — 128 re-executes `CatForce128` (built on demand, `-DCATFORCE_N=128`; rows are
  `unsigned __int128`, type `row_t` in LifeAPI.h). ~2.6-3.5x slower, identical suite results.
- Goal searches (any filter/match) also write `<output>_history.rle` in `[rule]History`
  (start pattern state 3, catalysts 1, goal cells 4).
- Layout files (`.rle`, multi-state, LifeHistory states; README "Layout files"): drawn roles
  1 pat, 2 catalyst-zone, 3 bait, 4 goal, 5 fixed-cat, 6 keep-empty, remappable with `legend`;
  `#C` keyword lines; `#CXRLE Pos` positions. `LayoutToLines` in CatForce.cpp turns a layout
  into keyword lines; the generator has a JS port (`layoutToLines`); both checked equal.
- `fixed-cat`, `catalyst-zone`, `no-catalyst-zone`, `keep-empty`, `torus-edge flag|drop|off`,
  filter ranges ending in `max`, `<<` RLE blocks, end-of-line `#` comments, `--show-input`.
- `split-level 2` — threads share second-level placements (helps many cores + few roots; slower on 2 cores).
- Catalyst options per `cat` line: `forbidden`, `required`, `antirequired`, `locus rle x y`,
  `mustinclude`, `transparent`, `check-recovery`, `sacrificial`, `period n`, `isblinker`
  (accept any phase of any oscillator).

Command line: `--check` (summary, no search), `--threads`, `--max-time`, `--max-results`, `--part`,
`--verbose`, `--compile-rule` (generates `rule_kernels.h` for the file's rule, rebuilds via `make`
or `$CATFORCE_MAKE`, restarts; keeps up to 8 rules), `--help`.

Engine changes (not visible in the input format):
- Compiled per-rule SIMD kernels (`RuleCompiler.h` → `rule_kernels.h`); table stepper fallback.
- Apple-silicon support via `sse2neon.h` and portable rotate macros.
- Mask caches `maskpack3-*` keyed by exact rule table + shapes + periods (safe to delete).
- Too-late placement pruning, duplicate removal by symmetry about any centre, envelope-sized output frames.
- Dynamic thread load balancing (threads claim the next top-level placement from a shared ticket).
- Lightcone pruning (§ above), contact filter, pause/resume.

## 4. Known facts and measurements worth keeping

- Profile of the user's 3-catalyst `3q4i-3` search: ~49% generation stepping (already SIMD and
  column-skipping), ~43% search bookkeeping/catalyst matching, ~6% candidate placement.
- Manual loci on 3q4i-3 (2-cat test): 190 s → 24 s but lost 4 of 8 results.
  `contact-filter yes` beat that (10.5 s, nothing lost).
- `max-active 100 → 20` gave ~40% time cut with identical results in an earlier test.
- Speed work (2026-10-17), results identical everywhere (suite, 2- and 3-catalyst checks):
  * Trajectories: each search node keeps its "no more catalysts" evolution; a child steps only
    the columns where it differs from its parent (narrowed every generation) and copies the rest
    in place. Falls back to full stepping when the difference covers >90% of live columns.
  * Too-late cutoff: a catalyst missing at a check generation g with g + ceil(S/p)*p > max-gen-1
    (S = stable-interval, p = period), or a missing bait with g + survive > max-gen-1, ends the
    branch (all descendants too).
  * Single-thread timings: HighLife max-active-20 search 10.2 s -> 5.7 s; 3q4i max-gen 100
    (max-active 100) 15.9 s -> ~12.7 s; 3q4i max-gen 35 15.9 s -> 6.3 s; 3-cat INT test 1.84 -> 0.78 s.
  * Profile before: ~60% of time was stepping leaf configurations (all catalysts placed).
  * Rejected: "board settled (period 1/2) with a catalyst missing" cutoff — fires often but the
    skipped generations were already nearly free under trajectories, and the per-gen hash cost
    more; skipping zero-init of temporaries / Move+Join scratch arrays — within noise.
  * Remaining profile (3q4i max-gen 100): ~37% full steps in leaves whose disturbance spread
    widely, 15% search bookkeeping, 12% column-limited steps, 6% each Convolve / MatchedBy.
- Things tried and rejected (no gain): passing masks by reference, PGO, break-after-success,
  dead-end pruning, incremental StepDirty, smaller lookup tables, memo caches,
  "remember single-catalyst outcomes" (no correct cheap version exists).
- Not done (expensive): full Lightcone-style approach signatures (5×5) per catalyst/orientation;
  a merge tool for `--part` / `_from<n>` outputs; resume that appends to the same output files.
- Diagnosis of the user's early `search.txt` (B3/S23-a5, D2| symmetry): results looked bad because
  (a) no goal/filter, (b) every `cat` used symmetry `.` (one orientation only), (c) num-catalyst
  means exactly N, (d) max-active 100. Corrected file in `user_files/search_fixed.txt`.

## 5. The generator (catforce_in_generator.html)

Built from three pieces in `generator/src/`:
- `a_head.html` — title, CSS (solid black, Arial, square corners, cool RGB-axis colours), markup.
- `c_ui.js` — all page logic (state, build/parse of .in, rendering, map, locus picker, drag regions).
- `old_generator_v8.html` — older page whose lines supply the shared helpers
  (cleanRle, decodeRle, checkRule, ruleTables, stepCells, classifyObject, splitComponents,
  clusterObjects, encodeRle, makeSetups).
- `build.py` — run `python3 build.py` inside `generator/src/` to rebuild
  `catforce_in_generator.html` (it also patches makeSetups to support the "Always present" role).

Generator features: six steps (rule, starting pattern with auto-split into active/bait/
always-present/leave-out objects, catalysts with stability check + orientation auto-detect +
locus picker on thumbnails, goals/filters, search size with drag-to-add regions on the map
(Shift = exclude), run/output options), warnings for under-oriented catalysts, long max-active,
no goal, etc. Draft is saved in browser storage. It was also published as a claude.ai artifact on
the old account; on a new account, publish `catforce_in_generator.html` again if wanted.

Generator now also has: torus size, spaceship-cutoff, torus-edge, contact-filter, lightcone,
unique-output, resume-from, split-level, fixed catalysts (+ import from an output-list line or a
result RLE), layout import/export, and an "Extra lines" box for keywords without a field.
`build.py` writes `generator/catforce_in_generator.html`; copy it to `CatForce/` too.
Page logic can be tested headless with jsdom (stub canvas getContext and clipboard).

## 6. Where things stood at the end

- Last user tasks: HighLife searches with an active 6/7-cell object plus several always-present
  objects. Working tested files: `user_files/highlife_fixed.txt` and `user_files/highlife2.txt`
  (need the current CatForce build because of the also-required fix).
- Generator: "Same search if closer than" and "Keep everything as one search" controls are now
  always shown after a paste; objects must be set to "Always present" individually.
- No open bugs known. Possible next steps the user discussed: test `contact-filter` on the full
  3-catalyst search; merge tool for split/resumed outputs; generator fields for the new options.

## 7. How to verify a build quickly

```
cd CatForce && make
./CatForce ../user_files/highlife2.txt --max-time 12     # ~1 s, 7 results
bash ../tests/suite.sh ./CatForce check                   # compare md5s with tests/base.txt
```
`tests/suite.sh` has hard-coded paths from the old machine: edit the `S=` line and the
`cp $S/ref/$f.in` line to point at this `tests/` folder before running.
