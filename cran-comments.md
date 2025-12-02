## Test environments
* local OS X install, R 4.5.1
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Response to CRAN Review
* **Non-standard files:** I have added 'docs', '_pkgdown.yml', and 'cran-comments.md' to .Rbuildignore to fix the top-level files note.
* **Misspelled words:** 'Borin' and 'Taglioni' are the surnames of the authors of the methodological papers cited. 'GVC' is the standard acronym for Global Value Chain, defined in the Description.
* **Resetting par():** Fixed all instances in vignettes where `par()` was changed without resetting. Now using `oldpar <- par(...)` at the start and `par(oldpar)` at the end of each visualization chunk in `bm2023-vs-bm2025.Rmd` and `bm2025-output-focus.Rmd`.

## Reverse dependencies
This is a new package, so there are no reverse dependencies.
