# GWmodel — vectorised, blocked, parallel GTWR

A fork of the CRAN package [GWmodel](https://cran.r-project.org/package=GWmodel),
built on release **2.4-1**. Same package name, same `library(GWmodel)`; the
internals of `gtwr()` and `st.dist()` are rewritten so that geographically and
temporally weighted regression can be fitted to a large panel.

Reference workload: a London LSOA-month panel of **34,013 rows**. On 15 cores it
fits in about 13 minutes at a 52 GB peak. Before this work it did not complete at
all above roughly 14,000 rows.

## What this changes relative to CRAN 2.4-1

| | CRAN 2.4-1 | Here |
|---|---|---|
| `st.dist` | six nested loops, cell by cell | vectorised, evaluated in column blocks |
| Unique-coordinate deduplication | absent | restored — 65 MB of lookups rather than an 8.6 GB expansion |
| `Q = (I−S)′(I−S)` | assembled in full, O(n³) | closed forms, O(n²) |
| Hat matrix `S` | materialised n×n | reductions accumulated per chunk |
| `gtwr()` parallelism | none; serial `for (i in 1:rp.n)` | `cores` argument, fork or PSOCK |
| `lamda` / `longlat` swap | live in `bw.gtwr.R` | fixed |

`bw.gtwr.R` in CRAN 2.4-1 calls `st.dist()` as `longlat = F, lamda = longlat`, so
`lamda` receives `FALSE` (zero) instead of its configured value. Since `lamda`
sets the spatial-versus-temporal mix, bandwidth selection weights distance
incorrectly whenever no distance matrix is supplied — with no error or warning.

## Verification

Every change is checked against reference fits rather than a tolerance. `st.dist`
output is bit-identical to the pre-optimisation implementation
(checksum `76668926238.3127899170` at n=4,000), and the reference GTWR fit
reproduces to 13 significant figures — the residue is upstream's rewrite of the
weighted least squares in `GWmodel.cpp`, which replaces an explicit `inv(xtwx)`
with a solve.

`bench/` holds the benchmark scripts and `tests/` the equivalence tests.

## Installing

Requires R ≥ 4.5 and a matching Rtools on Windows (Rtools45 covers 4.5 and 4.6);
the package contains C++ and will not install without a toolchain.

```r
install.packages("path/to/GWmodel", repos = NULL, type = "source")
```

---

## Authorship

This fork is a mix of my own work and AI-assisted work. The split:

**Mine, without AI**

- Identifying the GTWR performance bottleneck (2025).
- The original refactor from nested loops to vectorised arithmetic — the
  approach, and the first working implementation.
- Benchmarking and the equivalence tests (handover to AI at later point)

**AI-assisted (Claude Code), under my direction and review**

- Re-implementing that vectorisation in R.
- Additional fixes found in the process, most significantly **blocking** the
  distance construction and the related memory work: never materialising the hat
  matrix, replacing the O(n³) residual projector with closed forms, and giving
  each worker only its own columns of the distance matrix.
- The `lamda`/`longlat` correctness fix, and the Windows parallel path.


**Why the vectorisation was redone with AI rather than restored**

1. I no longer had the code for the original GTWR fix.
2. AI (Claude Opus 4.7) has become good enough that, with my oversight and my familiarity with this
   project, I was comfortable having it redo the work.

The bottleneck analysis and the decision to vectorise are mine and predate the AI
work; the blocking that made the full panel fit is not.
