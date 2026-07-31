## Test environments

* local R installation, R 4.6.1 (macOS, aarch64)
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

The one NOTE seen locally is:

```
* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  '__autograph_generated_filexu28yb1r.py' '__pycache__'
```

These files are written by Python's TensorFlow library (via its "autograph"
component) when the vignettes and tests execute model code. TensorFlow writes
them into the session temp directory and removes them from a Python `atexit`
handler, which is not always reached when Python is embedded in R.

This NOTE only appears on machines that have Python and TensorFlow installed.
Every test is guarded by `skip_if_not(check_tf_version())` and every vignette
chunk is guarded by `eval = check_tf_version("message")`, so on a machine
without TensorFlow no Python code runs at all and these files are never
created. The same NOTE was present for the `greta` 0.6.0 submission accepted in
July 2026.

## Submission notes

`greta.dynamics` was archived in September 2025 as a consequence of its
dependency, `greta`, being archived. `greta` 0.6.0 was accepted back onto CRAN
in July 2026, so this package can now be restored. This release requires
`greta` >= 0.6.0.

This version contains documentation fixes only; there are no user-facing
changes to package behaviour. The package has been tested against `greta` 0.6.0
and its full test suite passes.

## Method references

The methods implemented here are described in the `greta` paper already cited
in the DESCRIPTION, Golding (2019) <doi:10.21105/joss.01601>. There are no
additional publications describing this extension package specifically.

## revdepcheck results

There are no reverse dependencies.
