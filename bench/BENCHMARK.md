# GTWR vectorization benchmark

R 4.6.1 on Darwin 24.6.0 (arm64). Timings are `Sys.time()` deltas, median over reps. Old = master `R/gtwr.R` snapshot; new = current `R/gtwr.r` on `vectorize-gtwr`.

## What changed

- `st.dist()` `focus == 0` branches (both `rp.given` and symmetric): the two nested `for` loops became one vectorized matrix expression on a boolean-mask subset.
- `st.dist()` `focus > 0` branches: same idea, reduced to a single vector.
- `ti.distv()`: dropped the `c()`-growth loop; single vectorized `ti.dist()` call plus a mask for future observations.
- `ti.distm()`: kept the outer column loop (which is cheap) now that `ti.distv` is vector-native.
- `get.ts()`: replaced buggy `as.factor()` / `format()` round-trip with `sort(unique())` + `match()`. Also fixes a latent precision bug — see below.
- `get.uloat()`: replaced the O(n²) `which()` scan with a single `match()` on paste-keys.

## st.dist (called from `gtwr()` and `bw.gtwr()` per bandwidth candidate)

| Mode | n | old median | new median | speedup | equal | reps (old/new) |
|---|---:|---:|---:|---:|:---:|:---:|
| rp.given | 100 | 6.3 ms | 511 µs | 12.4× | yes | 20 / 20 |
| rp.given | 300 | 54.6 ms | 2.4 ms | 23.0× | yes | 10 / 10 |
| rp.given | 1000 | 578.0 ms | 23.3 ms | 24.8× | yes | 1 / 5 |
| rp.given | 2000 | 2.26 s | 72.5 ms | 31.1× | yes | 1 / 3 |
| symmetric | 100 | 3.7 ms | 563 µs | 6.6× | yes | 20 / 20 |
| symmetric | 300 | 31.1 ms | 3.0 ms | 10.4× | yes | 10 / 10 |
| symmetric | 1000 | 313.6 ms | 33.9 ms | 9.3× | yes | 1 / 5 |
| symmetric | 2000 | 1.30 s | 105.4 ms | 12.3× | yes | 1 / 3 |

## Helpers (called once per bandwidth trial and once per fit)

| Function | n | old median | new median | speedup | equal | reps (old/new) |
|---|---:|---:|---:|---:|:---:|:---:|
| get.ts n=5000 | — | 69.5 ms | 273 µs | 254.4× | **NO** | 5 / 5 |
| get.uloat n=5500 | — | 261.2 ms | 8.4 ms | 31.0× | yes | 5 / 5 |
| ti.distv n=5000 | — | 21.0 ms | 18 µs | 1160.2× | yes | 5 / 5 |

> `get.ts` shows `equal = **NO**` intentionally. The old implementation truncates full-precision numeric time stamps via `as.factor()` → `format()` (default `getOption("digits") = 7`), so `which(ts == i)` matches fewer entries than there are inputs and `index` comes back short. The new implementation returns a correct index of length `length(tv)`. See `tests/test_gtwr_equivalence.R` for the round-trip assertions.

## Notes on the shape of the speedups

- The `st.dist` speedup plateaus (~30× for `rp.given`, ~10–12× for symmetric at n ≥ 1000) once R interpreter overhead per iteration dominates the old path. Below n=100 setup costs on both sides pull the ratio down.
- The symmetric branch's ratio is lower than `rp.given` because the old code only iterated the upper triangle (~n²/2) while the new code touches the full n² before masking. One vectorized pass over n² still beats n²/2 R-level iterations by ~10×.
- Helper speedups (`ti.distv` ~800×, `get.uloat` ~30×, `get.ts` ~280×) reflect that the old code paid `c()`-growth allocation on every step. Under `gtwr()` these are called once per fit, so the wall-clock gain shows up mostly at large n or during bandwidth search where they run inside a golden-section loop.
- At n ≥ 1000 the old `st.dist` was timed with a single rep — one run takes 0.5–2.2 s and repeated timing added little signal; new-side medians are over 3–5 reps.
- Numbers exclude `s.dMat` / `t.dMat` construction, which is identical on both sides.

## Plots

Rendered by `bench/plot_results.R` after this script writes `bench/results.rds`.

- `bench/plots/st_dist_speedup.png` — speedup vs n for `st.dist`, one line per branch.
- `bench/plots/st_dist_wall_time.png` — log-log wall time, master vs vectorized, both branches.
- `bench/plots/helpers_speedup.png` — bar chart of helper speedups (`ti.distv`, `get.ts`, `get.uloat`).
- `bench/plots/parallel_speedup.png` — `.gtwr_dispatch` scaling from `bench/bench_parallel_synthetic.R` (no package install needed). The end-to-end `bench/bench_parallel.R` needs the compiled `GWmodel` package installed (which needs system GDAL/PROJ/GEOS); if you have those, it overwrites `parallel_results.rds` with real gtwr numbers.

## Correctness

See `tests/test_gtwr_equivalence.R` — 53 assertions across `st.dist` (all branches, master parity), helper round-trip properties, and the AICc guard. Runs in a couple of seconds via `Rscript tests/test_gtwr_equivalence.R`.

## Reproduce

```bash
Rscript bench/bench_st_dist.R      # this file — writes bench/results.rds
Rscript bench/bench_parallel.R     # optional — needs GWmodel + sp installed
Rscript bench/plot_results.R       # renders PNGs into bench/plots/
Rscript tests/test_gtwr_equivalence.R
```
