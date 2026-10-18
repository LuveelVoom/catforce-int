# CatForce
GoL Catalyst search utility based on LifeAPI library.

The main advantage of CatForce is that it doesn't make any assumptions
about the nature of the interaction. As long as the catalysts are back
in place in a given number of generations, they all could be destroyed
and reappear several times.

<!-- The torus centre is `(0, 0)` and left upper corner is `(-32, -->
<!-- -32)` and lower right corner is `(31,31)`. It has the same Y axis as -->
<!-- Golly (up is negative Y). -->

Usage
--
Run `make`, and then `./CatForce inputfile.in`. Currently only tested with `clang`.

Input File Format
--
See `examples/p83.in` etc. Some useful lists of catalysts are given in `catlists/`.

Parameters are separated by `" "` - i.e. space. Parentheses below
denote optional parameters.

| Line                  | Parameter                | Description                                                             |
|-----------------------|--------------------------|-------------------------------------------------------------------------|
| `max-gen`             | `n`                      | Overall maximum generation                                              |
| `start-gen`           | `n`                      | The min gen the _first_ encounter occurs                                |
| `last-gen`            | `n`                      | The max gen the _first_ encounter occurs                                |
| `num-catalyst`        | `n`                      | The number of catalysts to place (1 to 5)                               |
| `num-transparent`     | `n`                      | The number of catalysts marked `transparent` that may appear            |
| `stable-interval`     | `n`                      | Gens the catalysts must remain untouched to be considered stable        |
| `search-area`         | `x y w h`                | Search area; may be repeated, rectangles are combined (union)           |
| `exclude-area`        | `x y w h`                | Remove a rectangle from the search area; may be repeated                |
| `pat`                 | `rle`                    | The active pattern                                                      |
|                       | `(dx dy)`                | Offset applied to the active pattern                                    |
| `cat `                | `rle`                    | A catalyst                                                              |
|                       | `max-active`             | Number of generations in a row the catalyst may be missing              |
|                       | `dx dy`                  | Offset applied centre the catalyst, typically negative                  |
|                       | `symmetries-char`        | Character denoting the symmetry of the catalyst (see below)             |
|                       | `(forbidden rle x y)`    | Forbidden pattern around the catalyst (see below)                       |
|                       | `(required rle x y)`     | Cells of a catalyst that must stay ON in every generation               |
|                       | `(antirequired rle x y)` | Cells of a catalyst that must stay OFF in every generation              |
|                       | `(locus rle x y)`        | Cells of a catalyst that must interact first                            |
|                       | `(mustinclude)`          | Solutions must use at least one `mustinclude` catalyst                  |
|                       | `(transparent)`          | Marked as transparent for the purposes of `num-transparent`             |
|                       | `(check-recovery)`       | Always check the catalyst is recovered after exactly `max-active` gens  |
|                       | `(sacrificial)`          | Does not need to recover                                                |
|                       | `(period n)`             | Periodic catalyst                                                       |
|                       | `(isblinker)`            | Accept the catalyst in any phase (any oscillator, not just blinkers)    |
| `output`              | `filename`               | Output filename                                                         |
| `(and/or)filter`      | `genOrRange rle dx dy`   | Filter that must be matched for the solution to be accepted (see below) |
| `match`               | `genOrRange rle`         | Filter that checks whether the pattern occurs anywhere                  |
| `max-junk`            | `n`                      | For filters, the maximum non-matching non-catalyst population           |
| `full-report`         | `filename`               | Output filename for solutions ignoring all pattern filters              |
| `max-category-size`   | `n`                      | Maximum output row length before more solutions are dropped             |
| `output-gap`          | `n`                      | Empty cells between solutions in the output RLE (default 10)            |
| `keep-duplicates`     | `yes`                    | Disable duplicate removal (see below)                                   |
| `keep-interacting`    | `yes`                    | Keep solutions whose catalysts interact with each other (off by default) |
| `threads`             | `n`                      | Worker threads (default 0 = all cores); same results, just faster       |
| `part`                | `K/N`                    | Run only the K-th of N equal shares of the search (see below)           |
| `verbose`             | `yes`                    | Print every first-catalyst placement (the old console output)           |
| `max-time`            | `90`, `15m`, `2h`, `1d`   | Stop after this long; results so far are saved                          |
| `max-results`         | `n`                      | Stop after `n` solutions                                                |
| `unique-output`       | `filename`               | One square per distinct outcome, laid out in a square grid (default `<output>_unique.rle`) |
| `lightcone`           | `no`                     | Disable lightcone pruning (on by default) |
| `lightcone-tight`     | `yes`                    | Unsafe, faster: aim lightcone only at the mismatched cells at failure time (may lose results) |
| `resume-from`         | `n`                      | Skip the first n top-level placements (printed when a run is stopped; also saved to <input>.resume) |
| `contact-filter`      | `yes`                    | Skip placements where the new catalyst, with no further catalysts, never comes back within max-active gens (much faster; may lose results needing a rescue) |
| `split-level`         | `2`                      | Threads share out second-catalyst placements instead of first ones (better balance on many cores with few first placements; duplicates first-level work) |
| `output-list`         | `file`                   | Text list of every solution: catalyst shapes and positions              |
| `spaceship-cutoff`    | `n` or `off`             | Don't place catalysts for spaceships more than `n` cells from the rest of the reaction (default 10; see below) |
| `torus`               | `64` or `128`            | Torus size (default 64). `128` runs `CatForce128`, built automatically on first use (see below) |
| `fixed-cat`           | `rle x y max-active [cat options]` | A catalyst already in place (e.g. from an earlier result); see below |
| `catalyst-zone`       | `rle x y`                | Catalysts must lie entirely inside these cells; may be repeated (union) |
| `no-catalyst-zone`    | `rle x y`                | No catalyst cell may lie on these cells; may be repeated            |
| `keep-empty`          | `rle x y`                | These cells must stay dead in every generation; may be repeated     |
| `torus-edge`          | `flag`, `drop` or `off`  | Results whose reaction reaches the torus edge: marked (default), dropped, or not checked |
| `layout`              | `file.rle [x y]`         | Read a layout file (see **Layout files**) here                          |
| `legend`              | `4=goal:20-max ...`      | Change what layout states mean (for `layout` lines in this file)       |
| `include`             | `file`                   | Insert the lines of another file here (e.g. a catalyst list)            |
| `fit-in-width-height` | `w h`                    | Only allow solutions where all catalysts fit in a `w` by `h` rectangle  |
| `also-required`       | `rle x y`                | Require `rle` to be present in every generation                         |
| `bait`                | `rle x y [survive n] [absent n]` | Object in the starting pattern that must be disturbed, then reappear and survive |
| `symmetry`            | `symmetrycode`           | Global symmetry of the entire pattern (see below)                       |

**Catalyst Symmetry**: A character specifying what transformations are
applied to the catalyst:
- `|`: mirror by y.
- `-`: mirror by x.
- `+` or `@`: all 90 degree rotations, for D2- or D2/ invariant catalysts (eg boat).
- `x`: for 180 degree symmetrical catalysts.
- `/`: diagonal mirror.
- `*`: all 8 transformations.

**Filters**: Every potential solution will be checked against the
provided filters: only solutions that match the filters will be
kept. There are two variants: `andfilter` and `orfilter`. Every
`andfilter ...` must match and at least one `orfilter ...` must match
for the solution to succeed. A `match` filter pattern-matches the
solution for the specified pattern, and succeeds if the pattern occurs
anywhere (and otherwise behaves like an `andfilter`).

Filters that use a range of generations will only succeed after every
catalyst has interacted.

**Forbidden**: Checks for `rle` in the same location as the
catalyst. If `rle` is matched in any generation, the solution is
dropped. You can see `forbidden` as a "catalyst attached filter",
i.e. the filter location is based on the catalyst location. Intended
to exclude eaters/boat-bits and similar unwanted garbage. You may have
several forbidden patterns per catalyst. See the files in `catlists`
for examples.

**isblinker**: Normally a periodic catalyst must return in the same phase
it started in. With `isblinker` it counts as recovered if it matches *any*
of its phases exactly (cells alive, surrounding cells dead), so a blinker,
toad, beacon, etc. may come back rephased. Use it with `period n` (the
oscillator period); if `period` is omitted, the period is detected
automatically under the current rule.

**Compiled rule kernels (`--compile-rule`)**: for any rule other than B3/S23,
run once

    ./CatForce your.in --compile-rule

This writes `rule_kernels.h` with a stepping kernel generated for that rule (a
bit-parallel logic circuit: neighbour count, the rule's totalistic part, and
its exceptions such as 3q or 4i, minimised and run four columns at a time),
rebuilds the program with `make` (about 10-20 s), and then starts the search.
Searches in that rule are then about 2x faster (measured: B3/S23-q4i search
21.2 s -> 8.8 s, Cyclotronic 1.7 s -> 1.0 s; identical results). Up to 8 rules
are kept, so switching between a few rules needs no rebuild; a rule that is
not compiled in still works with the table stepper. The startup line says
which stepper is used, e.g. `(compiled kernel, 4.2x the table stepper here)`.
At startup each compiled kernel is checked against the table stepper on random
patterns (a mismatch disables it with a warning) and timed against it on your
machine; the faster one is used. Delete `rule_kernels.h` and run `make` to go
back to the table stepper only. If your build needs a different command than
`make` (e.g. `make CC=g++`), set it with the `CATFORCE_MAKE` environment
variable. The shipped `rule_kernels.h` already contains B3/S23-q4i and
Cyclotronic.

**Auto-split (generator page)**: paste a whole Golly pattern into the page's
"Auto-split" card. Each separate object is simulated alone under your rule;
still lifes and oscillators become `bait` lines, everything that moves or
evolves becomes the `pat`. Objects within a few cells of each other form one
setup, so a pattern containing several independent test cases becomes several
`.in` files (load one into the form, or copy/save all of them).

**Catalyst sanity check**: every solution's catalysts (including their mirror or
rotated images under `symmetry`) are evolved on their own for a full period and
must come back unchanged; otherwise the solution is rejected and counted in the
end summary. With a global `symmetry`, placements whose images would overlap or
touch the catalysts already placed are also skipped as soon as they are tried.
`keep-interacting yes` turns both off.

**Command line**: `./CatForce file.in [--check] [--show-input] [--threads N] [--max-time T]
[--max-results N] [--part K/N] [--verbose] [--compile-rule]`. `--check` parses the file, prints a
summary (rule, symmetry, pattern, number of catalysts and placements, search
positions, output files) and exits without searching; the other options
override the same keywords in the file.

**Late placements are skipped**: a solution needs every catalyst recovered for
`stable-interval` generations in a row before `max-gen`, so a catalyst placed
after generation `max-gen - stable-interval` can never be part of one. Those
placements are not tried (sacrificial catalysts excepted). Results are
unchanged; the saving is largest when `max-gen` is not much bigger than
`stable-interval` (measured on a B3/S23-q4i search: 2.7x at max-gen 30, 1.5x at
45, 1.2x at 60, with stable-interval 15).

**Hopeless branches are dropped**: once a catalyst is seen missing so late that
it can no longer stay recovered for `stable-interval` generations before
`max-gen` (or a bait can no longer survive long enough), that configuration and
everything that would be placed after it are abandoned. Results are unchanged;
it matters most with a large `max-active` (e.g. 100), where nothing else ends
these branches early.

**Console output**: by default the program prints one status line per second
(every 10 seconds after the first minute): results so far, elapsed time, and
how far the first-catalyst loop has got, e.g. `first catalyst at gen 12 of
1..60, 830 first placements done`. `--verbose` restores the old line for
every placement.

**Splitting a search** (`--part K/N`): runs only every N-th first-catalyst
placement, starting with the K-th, so N runs with `--part 1/N` ... `--part N/N`
together cover exactly the whole search (checked: the union equals a single
run, with no overlap). Use it to spread a search over several machines, or to
do a long search in pieces. Each part writes its own files, named with a
`_partKofN` suffix (e.g. `results_part2of4.rle`). Duplicate removal works within
a part, so mirror images can appear in two different parts.
Ctrl-C stops the search cleanly and saves everything found so far; press it a
second time to quit immediately.

**Threads**: all cores are used by default. The top-level placements are
split between the threads, so the set of solutions is the same as with one
thread (checked), only faster. Output order, and which of two mirror-image
solutions is kept, may vary between runs; `threads 1` makes runs fully
repeatable. Build with `-pthread` (the Makefile does).

**Output files**: results files start with `#C` comment lines (shown by Golly
under File > Show Comments) giving the number of solutions, the input file,
rule and settings, and whether the search had finished when the file was
saved. Files are written to a temporary name and then renamed, so stopping
the program during a save never leaves a half-written file.

**Include files**: `include catlists/The_Short_List_of_Catalysts` pulls in that
file's lines (relative to the including file, then the working directory;
nesting is allowed, Windows line endings are fine).

**Solution list**: with `output-list list.txt` each solution gets a line such as
`done by gen 26 | cat line 2 [<rle>] as <oriented rle> top-left (x,y) | ...`
so you can see where every catalyst sits without decoding the picture. At the
end the program prints how many solutions used each `cat` line.

**Mask cache**: the first run for a catalyst list writes a `maskpack3-*` file
in the working directory (collision masks between catalyst pairs). It is keyed
by the exact rule table, the exact catalyst shapes in order, their periods and
the options that skip pairs, so changing any of these recomputes it. Files are
written via a temporary name, so an interrupted run cannot leave a truncated
cache behind, and a damaged file is detected and recomputed. Create a
directory named `masks` to also cache per-pair masks (opt-in). Older
`maskpack-*`, `maskpack2-*` and `masks/maskraw*` files are ignored and can be
deleted; the oldest ones were not keyed by the rule and could give wrong
results.

**Input checking**: unknown keywords (typos) produce a warning instead of
being silently ignored; lines starting with `#`, `=` or `-` are comments.
Unknown `symmetry` names are an error (they used to mean "no symmetry"
silently; the old aliases `horizontal` and `rotate180evenboth` were affected
and now work). A catalyst's symmetry character must be one of `. | - \ / + x *`.
Numbers must be whole numbers, RLEs may only contain
`o b $` digits and `!` (no `x = ...` header, no spaces), and lines with
missing arguments are rejected with a message instead of being misread.

**Search area shapes**: the allowed catalyst positions are all `search-area`
rectangles combined, minus all `exclude-area` rectangles (with no `search-area`
line the default `-30 -30 60 60` is used). For example, to leave the top-right
quadrant empty:

    search-area -30 -30 60 60
    exclude-area 0 -32 32 32

The area limits where a catalyst's origin (the point set by its `dx dy`) may be
placed, not its whole extent. Coordinates: centre of the 64x64 torus (or 128x128 with `torus 128`) is 0,0;
x grows right, y grows down. The number of allowed positions is printed at
startup, and the search stops with a message if the area ends up empty.

**Bait**: `bait rle x y [survive n] [absent n]` adds an object to the
starting pattern (you do not need to repeat it in `pat`; position is given
like `pat`). A solution only counts if, together with all catalysts being
recovered, the bait is *missing* for at least `absent` generations (default 1,
so it must really be disturbed; use `absent 0` to also allow an untouched
bait) and is afterwards present again for `survive` consecutive generations
(default `stable-interval`). "Present" means its cells are alive and the
ring of cells around it is dead; oscillators are accepted in any phase (the
period is detected under the current rule). Under a global `symmetry`, the
bait is mirrored/rotated like the pattern and all images must reappear.
Several `bait` lines (up to 8) must all succeed. Catalysts may also be placed
to interact with the bait. At least one `cat` is still needed; the bait is
not a catalyst and is not counted in `num-catalyst`.

**Duplicates**: Solutions with an identical starting state are reported
only once, as are solutions that are mirror/rotated images of an
already-found one, when the active pattern (with any baits and
`also-required` cells) has that symmetry about *any* centre and the rule is
isotropic. The number removed is printed at the end. Add
`keep-duplicates yes` to keep everything.

**Writing input files**:
- `#` after a space starts a comment, also at the end of a line: `num-catalyst 2  # two`.
- `<<` in place of any RLE takes the RLE from the following lines, up to its `!`,
  exactly as Golly copies it (the `x = ...` header and `#` lines are skipped,
  line breaks are fine). The RLE keeps its own top-left corner, so offsets mean
  the same as on one line:

      pat << -3 -2
      x = 3, y = 4, rule = B3/S23
      o$3o$obo$
      2bo!

- A filter range may end in `max` (= `max-gen` − 1): `andfilter 20-max 2o$2o! 5 5`.
- `./CatForce file --show-input` prints the plain keyword lines a file amounts to
  (layouts and `<<` blocks expanded) and exits.

**Fixed catalysts (`fixed-cat`)**: `fixed-cat <rle> <x> <y> <max-active> [options]`
puts a catalyst in place before the search starts, with its top-left corner at
(x, y). It must recover like any other catalyst (the same `cat` options work:
`period`, `isblinker`, `required`, `forbidden`, ...), but it is not searched for
and does not count towards `num-catalyst`. Use it to search in stages: take a
result, fix its catalysts, and search for one more. Searched plus fixed
catalysts may be at most 8. The output list marks them `fixed-cat [...]`.
(Checked: fixing one catalyst of a 2-catalyst search and searching for 1 gives
exactly the 2-catalyst results that contain it.) The generator page can make
these lines from a pasted result.

**Torus edge (`torus-edge`)**: a result whose reaction (everything it touches
in its first `max-gen` + 14 generations) reaches the edge of the torus may have
wrapped around, which cannot happen on the plane. By default such results are
kept and marked `REACHES THE TORUS EDGE` in the `output-list` file, and the
count is printed at the end; `torus-edge drop` removes them, `torus-edge off`
skips the check. `torus 128` gives more room.

Layout files
---

A layout describes the geometry of a search as a picture: a multi-state RLE in
which each cell's state says what the cell is. Using LifeHistory states, you can
draw and edit it in Golly (rule `<your rule>History`, e.g. `B3/S23-q4iHistory`),
and the whole search (settings included) can live in one file that Golly opens
and CatForce runs directly:

    ./CatForce search.rle

Any file ending in `.rle`, or whose first line (after `#` lines) is an RLE
header, is read as a layout. A keyword file can also pull one in with
`layout file.rle [x y]`.

**What the states mean** (defaults; LifeHistory colours in brackets):

| State | Role | Becomes |
|---|---|---|
| 1 (green) | `pat`: the active pattern | one `pat` line |
| 2 (blue) | `zone`: catalysts must lie entirely inside these cells | `catalyst-zone` |
| 3 (white) | `bait`: each separate object is one bait | `bait` lines |
| 4 (red) | `goal`: each separate object must appear at exactly this spot | `andfilter 1-max ...` lines |
| 5 | `fixed`: each separate object is a catalyst already in place | `fixed-cat ... 20` lines |
| 6 | `keep-empty`: must stay dead in every generation | `keep-empty` |

With a zone and no `search-area` line, the whole torus is searched and the zone
alone limits where catalysts go; with both, both apply.

Odd states are alive at generation 0 and even states are dead, as in
LifeHistory, so the pattern, baits and fixed catalysts are drawn in live states
and the zone, goal and keep-empty cells in dead ones. "Separate objects" are
8-connected groups of cells of that state.

**Legend**: a `#C legend` line in the file (or a `legend` line in the keyword
file that has the `layout` line) changes the meaning of states, and gives roles
their options:

    #C legend 4=goal:30-max 5=fixed:15 3=bait:20:2 7=required 8=exclude 9=ignore

Roles: `pat`, `zone`, `exclude` (no catalyst cell may lie here), `bait[:survive[:absent]]`,
`required` (always alive, `also-required`), `fixed[:max-active[:cat options]]`,
`goal[:generations]` (default `1-max`), `keep-empty`, `ignore`. States up to 24
(letters A to X) can be given roles; a state that appears in the drawing with no
role is an error, so a stray state is never silently ignored.

**Settings**: `#C` lines whose first word is a CatForce keyword are read as
keyword lines (`#C max-gen 100`, `#C cat 2o$2o! 20 -1 -1 .`, `#C include catlist`);
other `#C` lines are ordinary comments. Plain keyword lines may also follow the
RLE after its `!`. Keyword lines are read after the drawing, so they can add to
it (more `bait` or `cat` lines), and a `rule` line overrides the header.

**Rule**: taken from the header (`rule = B3/S23-q4iHistory` gives `B3/S23-q4i`;
a trailing `History` or `Super` is removed), unless a `rule` line is given.

**Position**: the layout's top-left cell goes to the `#CXRLE Pos=x,y` that Golly
writes (so coordinates match Golly's), else to the `x y` of the `layout` line,
else the layout is centred on (0, 0). The first line of `--show-input` says where
it went. It must fit on the torus (use `torus 128` for big layouts).

`examples/layout_example.rle` is a small complete example, and the generator
page can save the current search as a layout.

**Escaped spaceships (`spaceship-cutoff`)**: a spaceship that has left the
reaction (more than `n` cells from every other non-catalyst cell, with the
reaction or its debris still present) no longer triggers catalyst placements,
and the positions it passes are blocked for later catalysts (one placed there
would have been hit by it). This removes the "extra conduits" where a later
catalyst just eats or deflects an escaping glider far from the reaction, or a
spaceship that wrapped around the torus. A cluster counts as a spaceship if it
comes back moved within 32 generations when run alone. The default is 10;
`spaceship-cutoff off` restores the old behaviour. Note that a spaceship in the
starting pattern that is already more than `n` cells from the rest of the
pattern is treated the same way.

**History output (goal searches)**: when the file has any `filter`,
`andfilter`, `orfilter` or `match` line, a third file `<output>_history.rle` is
written in the rule's `History` variant (e.g. `B3/S23-q4iHistory`, open it in
Golly 4.0 or later). It has the same layout as the main output; the starting
pattern is in state 3 (white), so it leaves red cells where it dies, the
catalysts are state 1 (green), and the cells where the goal pattern was found
are state 4 (red).

**Torus size (`torus 128`)**: the default torus is 64x64 (coordinates -32..31).
With `torus 128` the search runs on a 128x128 torus (coordinates -64..63), so
search areas can be larger and spaceships take twice as long to wrap around.
`CatForce` then starts `CatForce128` from the same folder, building it with
`make CatForce128` the first time (afterwards plain `make` and
`--compile-rule` keep it up to date). Results are the same as on the 64 torus
when nothing reaches the edge (checked on the test suite). It is slower:
measured 2.6-3.5x per search on the test files (8x with `match` filters), and
building the mask cache the first time for a catalyst list takes about 16x
as long.

**Output layout**: each row of solutions is drawn in a shared frame sized to
everything those solutions touch while they run (not just their starting
cells), with `output-gap` empty cells between frames. So the output can be
run in Golly without neighbouring solutions colliding (spaceships that escape
the reaction can still reach a neighbour eventually).

**Check Recovery**: Catalysts with the `check-recovery` attribute are
immediately tested for whether they recover in `max-active`
generations, without the support of any further catalysts. This is
useful for catalysts with a bait still life that has a long recovery
time, like the hive-pushes or the loaf-spin catalysts.

<!-- Combining Results -->
<!-- --- -->

<!-- `combine-results yes [<survive-0> <survive-1> ...]` -->

<!-- If this feature is enabled the search will at first ignore all filters -->
<!-- and survival inputs, and will search all the possible catalysts. Then -->
<!-- it will try to combine all the found catalysts in all possible -->
<!-- combinations, and only then will filter by `survive-i` and apply the -->
<!-- filters to exclude them from the final report. -->

<!-- This feature will generate report as follows: -->

<!-- - `output.rle` - all the possible catalysts. -->
<!-- - `output.rle_Combined*.rle` - will generate all combined reports. -->
<!-- - `output.rle_Final.rle` - the final report. **This is the main output.** -->

<!-- Optional survival filter per "iteration" are added. Combine works as -->
<!-- follows: each time it start from the initial search results (combine -->
<!-- by default uses survive count = 1), and tries to add catalyst from -->
<!-- those results. Sometimes one could get explosion, if the interaction -->
<!-- is very potent. So filter is added to limit the combine, by surviving -->
<!-- count (if something doesn't survive with two catalyst for 5 -->
<!-- iterations, it's probably junk - so CatForce will filter it on the -->
<!-- second combine iteration and not in the end). -->

<!-- This allows faster and more efficient combine operation with very -->
<!-- potent conduits which otherwise would overflow the system, with many -->
<!-- useless catalysts. -->

<!-- **NOTE** Recommended for use only for `num-catalyst` = 1 or 2 -->

<!-- **NOTE** See 4.in file for example.  -->

<!-- **NOTE** CatForce will use the last `survive-i` as the default from -->
<!-- that point on. If you don't enter any numbers it will use survival -->
<!-- count 1, and will filter only when finish all possible combinations. -->

Choosing the Rule (INT rules)
---

Put a `rule` line anywhere in the `.in` file; no code changes or
recompiling are needed to switch rules:

    rule B3/S23                       (Conway's Life; also the default if omitted)
    rule B3ainj4cw5cek6n/S2-in3-ky4enyjtw5injr6cen7c   (Cyclotronic)
    rule B2ei3aeikqr/S2-a3            (any Hensel isotropic non-totalistic rule)
    rule MAPARYXfhZofugWaH7oaIDogBZofuhogOiAaIDogIAAgAAWaH7oaIDogGiA6ICAAIAAaIDogIAAgACAAIAAAAAAAA

Accepted formats are Life-like / Hensel strings (`B.../S...`, `S/B` numeric
form, upper or lower case, letters and `-` negation as in Golly) and Golly
`MAP` strings. Invalid rules stop the program with an explanatory message,
and the parsed rule is echoed at startup. Rules containing B0 are not
supported. B3/S23 uses the original fast stepper; every other rule uses a
table-driven stepper, which updates 4 cells per lookup in a 256 KB table and
visits only the parts of the board next to live cells; per generation it is
roughly half the speed of the Conway stepper. `--compile-rule` (see below)
builds a kernel for your rule that is about as fast as the Conway one. Non-isotropic `MAP` rules are supported (checked
against Golly's MAP definition); for them, rotated or mirrored catalyst
variants may behave differently and duplicate removal by symmetry is off.
Output files carry the rule as given (an isotropic `MAP` rule is written as
its Hensel equivalent). Build with `make` (or `make CC=g++` if you do not
have clang++).

Symmetric Searches
---


`symmetry s` automatically applies the symmetry `s` to the active
pattern and catalysts. The options for `s` are:
- D2 symmetries: `D2|` (reflect across x = 0), `D2-` (y = 0), `D2\` (y=x) and `D2/` (y=-x)
- C2 symmetries: `C2_1` (bounding box odd by odd), `C2_2` (odd by even), `C2_4` (even by even)
- C4 symmetries: `C4_1`, (odd by odd) `C4_4` (even by even)
- D4 symmetries: `D4_+1`, `D4_+2` (odd by even),`D4_+4`, `D4_x1`, and `D4_x4`
- D8 symmetries: `D8_1` and `D8_4`

The options D2 follow
[LLS](https://gitlab.com/OscarCunningham/logic-life-search), while the
rest of the symmetries follow
[Catagolue](https://catagolue.hatsya.com/census). [LifeWiki](https://conwaylife.com/wiki/Static_symmetry)
has a more in-depth explanation, though the notation there is slightly
different. (They choose `C2_2` to be even by odd not odd by even).

There are LLS-inspired alternatives for the other groups:
- D2 symmetries: as above
- C2 symmetries: `C2 C2even C2|even C2-even`
- C4 symmetries: `C4 C4even`
- D4 symmetries: `D4+ D4+|even D4+-even D4+|even D4x D4xeven`
- D8 symmetries: `D8 D8even`

The default is "odd" transformations, i.e. those that fix the cell at
(0,0); `verticaleven` and `horizontaleven` (referring to bounding box
dimensions) are recognized as equivalent to `-even` and `|even`,
respectively.
