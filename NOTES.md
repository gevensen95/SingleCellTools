# Open investigations

## RCTD: `RunRCTD()` fails with `assignInNamespace` lock error (unresolved)

**Symptom:**

```r
RunRCTD(visium, assay_query = 'SCT', assay_ref = 'SCT', reference = sennet,
        celltype_col = "celltype", mode = "full",
        max_cells_per_ref_celltype = 5000, n_cores = 10)
#> Error in utils::assignInNamespace("detectCores", ...) :
#>   locked binding of 'detectCores' cannot be changed
```

**Where it comes from:** `R/RunRCTD.R`'s `.run_rctd_one()` has a defensive patch
that monkey-patches `parallel::detectCores()` via `assignInNamespace()`. It was
added because `spacexr`'s internals call `parallel::detectCores()` unguarded
and crash on `NA`. In this environment, the namespace binding can't be
unlocked, so the patch itself now errors before RCTD even runs.

**Root cause, found via two diagnostics:**

```r
> sessionInfo()$R.version.string
Error in system(paste(which, shQuote(names[i])), intern = TRUE, ignore.stderr = TRUE) :
  cannot popen '/bin/which 'timedatectl' 2>/dev/null', probable reason 'Cannot allocate memory'
> parallel::detectCores()
[1] NA
```

The R process can't `fork()`/`popen()` at all in this environment — this is
not a narrow `detectCores()` quirk, it's a fork/subprocess restriction at the
OS or container level. That means `spacexr`'s internal `mclapply`-based
multicore execution will likely also fail regardless of any `detectCores()`
patch — `n_cores > 1` may not be viable here at all.

**Two options, not yet decided between:**

1. Drop the `assignInNamespace` patch entirely and check whether
   `detectCores()` behaves differently in a fresh R session (maybe the lock
   is session-state, not environment-wide).
2. Force serial execution (`n_cores = 1`) given the environment can't fork,
   and/or shim the underlying OS probe (e.g. stub out whatever calls
   `/bin/which timedatectl`) rather than touching R's locked namespace.

Next step: reproduce in a clean R session to see if (1) is viable before
committing to (2).
