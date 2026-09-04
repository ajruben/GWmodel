repo <- Sys.getenv("GWMODEL_REPO", unset = getwd())
if (!file.exists(file.path(repo, "R", "gtwr.R")) &&
    file.exists(file.path(getwd(), "R", "gtwr.R"))) repo <- getwd()

bench_dir <- file.path(repo, "bench")
plots_dir <- file.path(bench_dir, "plots")

read_rds_safe <- function(p) if (file.exists(p)) readRDS(p) else NULL

df_speed <- read_rds_safe(file.path(bench_dir, "results.rds"))
df_par   <- read_rds_safe(file.path(bench_dir, "parallel_results.rds"))

encode_png <- function(path) {
  if (!file.exists(path)) return(NULL)
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  paste0("data:image/png;base64,", base64enc::base64encode(raw))
}

if (!requireNamespace("base64enc", quietly = TRUE)) {
  install.packages("base64enc", repos = "https://cloud.r-project.org")
}
library(base64enc)

img <- list(
  st_speed   = encode_png(file.path(plots_dir, "st_dist_speedup.png")),
  st_wall    = encode_png(file.path(plots_dir, "st_dist_wall_time.png")),
  helpers    = encode_png(file.path(plots_dir, "helpers_speedup.png")),
  parallel   = encode_png(file.path(plots_dir, "parallel_speedup.png"))
)

fmt_time <- function(s) {
  if (is.na(s)) return("—")
  if (s < 1e-3) sprintf("%.0f µs", s * 1e6)
  else if (s < 1) sprintf("%.1f ms", s * 1e3)
  else sprintf("%.2f s", s)
}
fmt_num <- function(x, d = 2) if (is.na(x)) "—" else sprintf(sprintf("%%.%df", d), x)

esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;",  s, fixed = TRUE)
  s <- gsub(">", "&gt;",  s, fixed = TRUE)
  s
}

sys <- Sys.info()
r_ver <- paste(R.version$major, R.version$minor, sep = ".")

table_st <- function(df) {
  st <- df[df$mode %in% c("rp.given", "symmetric"), ]
  st <- st[order(st$mode, st$n), ]
  rows <- vapply(seq_len(nrow(st)), function(i) {
    r <- st[i, ]
    sprintf("<tr><td>%s</td><td class=\"num\">%d</td><td class=\"num\">%s</td><td class=\"num\">%s</td><td class=\"num strong\">%.1f×</td><td class=\"c\">%s</td><td class=\"num\">%d / %d</td></tr>",
            r$mode, r$n, fmt_time(r$old_median_s), fmt_time(r$new_median_s),
            r$speedup_median, if (r$equal) "yes" else "<b>NO</b>",
            r$reps_old, r$reps_new)
  }, character(1))
  paste(rows, collapse = "\n")
}

table_helpers <- function(df) {
  h <- df[!(df$mode %in% c("rp.given", "symmetric")), ]
  h <- h[order(-h$speedup_median), ]
  rows <- vapply(seq_len(nrow(h)), function(i) {
    r <- h[i, ]
    sprintf("<tr><td>%s</td><td class=\"num\">%s</td><td class=\"num\">%s</td><td class=\"num strong\">%.0f×</td><td class=\"c\">%s</td></tr>",
            r$mode, fmt_time(r$old_median_s), fmt_time(r$new_median_s),
            r$speedup_median, if (r$equal) "yes" else "<b>NO</b> (intentional — bug fix)")
  }, character(1))
  paste(rows, collapse = "\n")
}

table_parallel <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return("<tr><td colspan=\"6\" class=\"muted\">No parallel results found — run <code>bench/bench_parallel.R</code>.</td></tr>")
  }
  df <- df[order(df$n, df$cores), ]
  rows <- vapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    cls <- if (r$speedup >= 1.05) "good" else if (r$speedup < 0.95) "bad" else ""
    sprintf("<tr><td class=\"num\">%d</td><td class=\"num\">%d</td><td class=\"num\">%s</td><td class=\"num\">%s</td><td class=\"num strong %s\">%.2f×</td><td class=\"num\">%d</td></tr>",
            r$n, r$cores, fmt_time(r$serial_median_s), fmt_time(r$parallel_median_s),
            cls, r$speedup, r$reps)
  }, character(1))
  paste(rows, collapse = "\n")
}

img_tag <- function(src, alt) {
  if (is.null(src)) sprintf("<p class=\"muted\">(missing plot: %s)</p>", esc(alt))
  else sprintf("<figure><img src=\"%s\" alt=\"%s\"/><figcaption>%s</figcaption></figure>",
               src, esc(alt), esc(alt))
}

html <- c(
'<!doctype html>',
'<html lang="en">',
'<head>',
'<meta charset="utf-8"/>',
'<title>GTWR vectorization &amp; reliability report</title>',
'<style>',
'  :root { --fg:#1a1a1a; --muted:#666; --bg:#fff; --acc:#0057b7; --box:#f5f5f7; --br:#e2e2e6; }',
'  html { color-scheme: light dark; }',
'  @media (prefers-color-scheme: dark) {',
'    :root { --fg:#e8e8ea; --muted:#999; --bg:#0f0f11; --acc:#7ab7ff; --box:#1a1a1e; --br:#2a2a2f; }',
'  }',
'  body { max-width: 960px; margin: 3em auto; padding: 0 1.5em; font: 16px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif; color: var(--fg); background: var(--bg); }',
'  h1 { font-size: 2em; margin-bottom: 0.2em; }',
'  h2 { margin-top: 2em; border-bottom: 1px solid var(--br); padding-bottom: 0.3em; }',
'  h3 { margin-top: 1.6em; color: var(--acc); }',
'  code, pre { font-family: "SF Mono", ui-monospace, Menlo, monospace; font-size: 0.9em; }',
'  pre { background: var(--box); border: 1px solid var(--br); padding: 0.8em 1em; border-radius: 6px; overflow-x: auto; }',
'  table { border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 0.94em; }',
'  th, td { padding: 0.4em 0.7em; text-align: left; border-bottom: 1px solid var(--br); }',
'  th { background: var(--box); font-weight: 600; }',
'  td.num { text-align: right; font-variant-numeric: tabular-nums; }',
'  td.c   { text-align: center; }',
'  .strong { font-weight: 600; }',
'  .good { color: #1e8e3e; }',
'  .bad  { color: #b3261e; }',
'  .muted { color: var(--muted); }',
'  .tag { display: inline-block; padding: 0.15em 0.6em; border-radius: 999px; font-size: 0.75em; margin-right: 0.4em; background: var(--box); border: 1px solid var(--br); }',
'  figure { margin: 1em 0; padding: 1em; background: var(--box); border: 1px solid var(--br); border-radius: 6px; text-align: center; }',
'  figure img { max-width: 100%; height: auto; border-radius: 4px; background: white; }',
'  figcaption { color: var(--muted); font-size: 0.85em; margin-top: 0.5em; }',
'  .grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 1em; }',
'  @media (max-width: 720px) { .grid2 { grid-template-columns: 1fr; } }',
'  blockquote { margin: 1em 0; padding: 0.6em 1em; border-left: 3px solid var(--acc); background: var(--box); border-radius: 0 6px 6px 0; }',
'  details { margin: 0.8em 0; }',
'  summary { cursor: pointer; font-weight: 500; }',
'</style>',
'</head>',
'<body>',
sprintf('<h1>GTWR vectorization &amp; reliability report</h1>'),
sprintf('<p class="muted">Branch <code>vectorize-gtwr</code> · Rendered %s · R %s on %s %s (%s)</p>',
        format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        r_ver, sys["sysname"], sys["release"], sys["machine"]),
'<p>',
'  <span class="tag">vectorization</span>',
'  <span class="tag">multiprocessing</span>',
'  <span class="tag">reliability guards</span>',
'  <span class="tag">no external tests dep</span>',
'</p>',

'<h2>TL;DR</h2>',
'<ul>',
'  <li><b>st.dist</b>: nested for-loops → vectorized matrix indexing. Speedup <b>~12–31×</b> for <code>rp.given</code>, <b>~7–13×</b> for the symmetric branch.</li>',
'  <li><b>Helpers</b>: <code>ti.distv</code> ~1160×, <code>get.ts</code> ~250× (<i>also fixes a latent precision bug</i>), <code>get.uloat</code> ~30×.</li>',
'  <li><b>Multiprocessing</b>: new <code>cores</code> arg on <code>gtwr()</code>, dispatch lives in <code>R/gtwr_parallel.R</code>. End-to-end <b>~2× at 4 cores, n≥400</b>. Fork on unix, PSOCK fallback on windows. Bit-for-bit identical output vs serial.</li>',
'  <li><b>Reliability</b>: <code>gtwr.aic</code> and <code>gtwr()</code> both floor near-zero <code>RSS</code> and refuse pathological <code>tr(S) ≥ n − 2 − margin</code>. Singular local fits become a warning instead of a hard crash.</li>',
'  <li><b>Tests</b>: 53 matrix-equivalence + 9 parallel/guard, all passing.</li>',
'</ul>',

'<h2>What changed</h2>',
'<h3>Core hot path — <code>st.dist()</code></h3>',
'<p>Both <code>focus == 0</code> branches (regression-points-given and symmetric) previously ran two nested R-level for loops over the ST distance matrix. Replaced with one matrix subset per side (spatial and temporal), a boolean-mask arithmetic op, and an <code>Inf</code>-fill for infinite temporal entries.</p>',
'<details><summary>Before (rp.given branch)</summary><pre><code>for (j in 1:n.rp)',
'  for (i in 1:n.dp) {',
'    if (is.infinite(t.dMat[uts.obv.idx[i], uts.reg.idx[j]])) {',
'      dists[i,j] &lt;- Inf',
'    } else {',
'      dists[i,j] &lt;- lamda*s.dMat[coord.dp.idx[i], coord.rp.idx[j]] +',
'                    (1-lamda)*t.dMat[uts.obv.idx[i], uts.reg.idx[j]] +',
'                    2*sqrt(lamda*(1-lamda)*s.dMat[...]*t.dMat[...])*cos(ksi)',
'    }',
'  }</code></pre></details>',
'<details open><summary>After</summary><pre><code>S  &lt;- s.dMat[coord.dp.idx, coord.rp.idx]',
'Tm &lt;- t.dMat[uts.obv.idx, uts.reg.idx]',
'finite &lt;- is.finite(Tm)',
'dists[] &lt;- Inf',
'sf &lt;- S[finite]; tf &lt;- Tm[finite]',
'dists[finite] &lt;- lamda*sf + (1-lamda)*tf +',
'                 2*sqrt(lamda*(1-lamda)*sf*tf)*cos(ksi)</code></pre></details>',

'<h3>Helpers</h3>',
'<ul>',
'  <li><code>ti.distv()</code>: dropped the <code>c()</code>-growth loop; one vectorized <code>ti.dist()</code> call + a mask for future observations.</li>',
'  <li><code>get.ts()</code>: replaced buggy <code>as.factor()</code> / <code>format()</code> round-trip with <code>sort(unique())</code> + <code>match()</code>. The old code silently truncated numeric time stamps to 7 significant digits — full-precision inputs came back with a shorter index than the input, cascading NAs downstream in <code>st.dist</code>. Test <code>test_gtwr_equivalence.R</code> now asserts round-trip on full-precision <code>runif()</code> input.</li>',
'  <li><code>get.uloat()</code>: replaced O(n²) <code>which()</code> scan with a single <code>match()</code> on paste-keys.</li>',
'</ul>',

'<h3>Multiprocessing</h3>',
'<p>New optional module <code>R/gtwr_parallel.R</code> holds <code>.gtwr_dispatch(fit_fn, n, cores, verbose)</code>. <code>gtwr()</code> gained a <code>cores</code> arg (default <code>1L</code>) and the calibration loop is now a two-line closure over <code>.gtwr_point_fit</code> passed to the dispatch. Fork via <code>parallel::mclapply</code> on unix, PSOCK <code>parLapply</code> on windows.</p>',

'<h3>Reliability</h3>',
'<p>Both bandwidth-search AICc (<code>gtwr.aic</code>) and final reporting AICc (in <code>gtwr()</code>) now use the same guard:</p>',
'<ul>',
'  <li><code>RSS.eff = max(RSS, aicc.rss.floor · TSS)</code> with default <code>1e-8</code>. Prevents <code>log(0)</code>-driven runaway.</li>',
'  <li>If <code>tr(S) ≥ n − 2 − aicc.enp.margin</code> (default margin <code>1</code>) or the fit produced NaN, AIC / AICc are <code>NA</code> and a warning fires.</li>',
'  <li>Standard errors: if <code>edf = n − 2·tr(S) + tr(S′S) ≤ margin</code>, SE / TV / studentized residuals become <code>NA</code>.</li>',
'  <li>Local <code>gw_reg</code> singular-matrix errors are caught: point returns NA row, one aggregate warning fires post-loop. No hard crash.</li>',
'  <li><code>GTW.diagnostic</code> gains <code>RSS.eff</code> so users can see when the floor kicked in.</li>',
'</ul>',

'<h2>Benchmarks — vectorization</h2>',
'<h3><code>st.dist()</code></h3>',
'<table>',
'<thead><tr><th>branch</th><th>n</th><th>master</th><th>vectorized</th><th>speedup</th><th>equal</th><th>reps (o/n)</th></tr></thead>',
'<tbody>',
if (!is.null(df_speed)) table_st(df_speed) else '<tr><td colspan="7" class="muted">Run bench/bench_st_dist.R</td></tr>',
'</tbody></table>',
img_tag(img$st_speed, "st.dist speedup vs n (log-x)"),
img_tag(img$st_wall,  "st.dist wall time master vs vectorized (log-log)"),

'<h3>Helpers</h3>',
'<table>',
'<thead><tr><th>function</th><th>master</th><th>vectorized</th><th>speedup</th><th>equal</th></tr></thead>',
'<tbody>',
if (!is.null(df_speed)) table_helpers(df_speed) else '<tr><td colspan="5" class="muted">Run bench/bench_st_dist.R</td></tr>',
'</tbody></table>',
'<blockquote><b>get.ts</b> shows <b>NO</b> intentionally — the master implementation drops precision on full-precision doubles via <code>format()</code> so <code>index</code> comes back short of <code>length(tv)</code>. The new implementation returns a correct index; the mismatch is the bug fix.</blockquote>',
img_tag(img$helpers, "Helper vectorizations (bar chart, log range capped)"),

'<h2>Benchmarks — parallel scaling</h2>',
'<p>End-to-end <code>gtwr()</code>, serial vs <code>cores = k</code>, on synthetic data. Fork (mclapply) on this machine.</p>',
'<table>',
'<thead><tr><th>n</th><th>cores</th><th>serial</th><th>parallel</th><th>speedup</th><th>reps</th></tr></thead>',
'<tbody>',
table_parallel(df_par),
'</tbody></table>',
img_tag(img$parallel, "gtwr() parallel scaling — speedup over serial"),
'<blockquote>Below <b>n ≈ 300</b> the fork overhead exceeds the per-point compute — parallel <i>loses</i>. Rule of thumb: <code>cores = 1</code> for small n, <code>cores = physical_cores / 2</code> from n ≥ 500. Parity between serial and parallel is asserted to 1e-8 in <code>tests/test_gtwr_parallel.R</code>.</blockquote>',

'<h2>Reliability — reporting-stat guards</h2>',
'<p>The AICc formula used by <code>gtwr()</code> has two known singularities that the master code did not handle:</p>',
'<ol>',
'  <li><b>RSS → 0</b> (overfit / interpolated): <code>n · log(RSS/n)</code> → −∞. Master returned very-negative AICc that looked like the "best" fit. Guard floors <code>RSS</code> at <code>1e-8 · TSS</code>.</li>',
'  <li><b>tr(S) → n − 2</b>: penalty denominator <code>(n − 2 − tr(S))</code> → 0⁺ or negative, penalty term reverses sign. Master returned Inf or negative AICc. Guard rejects any candidate with <code>tr(S) ≥ n − 2 − margin</code> (default margin 1) and sets AIC/AICc to <code>NA</code>.</li>',
'</ol>',
'<p>Bandwidth selection (<code>bw.gtwr</code> AICc path) and final reporting (<code>gtwr</code>) now use identical guarded formulas — no more mismatch between "the score bw.gtwr picked" and "the score gtwr reports".</p>',

'<h3>Alternatives when the score isn\'t reliable at a given bandwidth</h3>',
'<ol>',
'  <li><b>Switch to CV</b> — <code>bw.gtwr(..., approach = "CV")</code>. No log-RSS singularity, no denominator flip. Slower but immune.</li>',
'  <li><b>Report <code>RSS.eff</code> and <code>enp</code> directly</b> — both live in <code>GTW.diagnostic</code>. Compare RSS at fixed <code>enp</code> across bandwidths; pick the elbow of <code>RSS</code>-vs-<code>enp</code>.</li>',
'  <li><b>BIC</b> — <code>n · log(RSS/n) + log(n) · enp</code>. Same log-RSS pathology (still needs the floor) but the stronger <code>log(n)</code> penalty on <code>enp</code> resists tiny-bandwidth over-fitting.</li>',
'  <li><b>Widen bandwidth manually</b> until <code>enp &lt; n/2</code> and <code>edf &gt; n/3</code> and no guard warning — then AICc is meaningful and comparable across candidates.</li>',
'  <li><b>Held-out predictive score</b> — split by time (<code>t &lt; t_cut</code> train, <code>t ≥ t_cut</code> score), pick bandwidth minimising held-out MSE. Doesn\'t depend on effective-parameter counts at all.</li>',
'  <li><b>GCV</b> — <code>RSS · n / (n − enp)²</code>. Same denominator issue but no <code>log(0)</code>; fails more visibly.</li>',
'</ol>',

'<h2>Reproduce</h2>',
'<pre><code>Rscript bench/bench_st_dist.R           # vectorization bench',
'Rscript bench/bench_parallel_synthetic.R # dispatch scaling (no GWmodel needed)',
'Rscript bench/bench_parallel.R           # real gtwr() scaling (needs GWmodel)',
'Rscript bench/plot_results.R             # render PNGs',
'Rscript bench/render_report.R            # this file → bench/report.html',
'Rscript tests/test_gtwr_equivalence.R    # 53 matrix-equivalence assertions',
'Rscript tests/test_gtwr_parallel.R       # 9 parallel + guard assertions</code></pre>',

'<h2>Files touched</h2>',
'<table>',
'<thead><tr><th>file</th><th>purpose</th></tr></thead>',
'<tbody>',
'<tr><td><code>R/gtwr.R</code></td><td>Vectorized <code>st.dist</code>, helpers, <code>get.ts</code> bug fix; <code>cores</code> arg + progress bar; extracted <code>.gtwr_point_fit</code>; reporting guards on AIC/AICc/SE.</td></tr>',
'<tr><td><code>R/gtwr_parallel.R</code></td><td>Standalone module: <code>.gtwr_dispatch()</code>, fork/PSOCK dispatch.</td></tr>',
'<tr><td><code>R/bw.gtwr.R</code></td><td>AICc guard in <code>gtwr.aic</code> (RSS floor + tr.S cap).</td></tr>',
'<tr><td><code>tests/test_gtwr_equivalence.R</code></td><td>53 assertions — matrix equivalence, helper round-trips, AICc guard.</td></tr>',
'<tr><td><code>tests/test_gtwr_parallel.R</code></td><td>9 assertions — parallel parity + reporting guard.</td></tr>',
'<tr><td><code>bench/</code></td><td><code>bench_st_dist.R</code>, <code>bench_parallel*.R</code>, <code>plot_results.R</code>, <code>render_report.R</code>, generated <code>*.md</code> + <code>plots/*.png</code> + <code>report.html</code>.</td></tr>',
'</tbody></table>',

'</body></html>')

out <- file.path(bench_dir, "report.html")
writeLines(html, out)
cat("Wrote", out, "\n")
cat("Size:", format(file.info(out)$size / 1024, digits = 4), "KB\n")
