## Submission

This is an update of SynthesizerPlus (previous submission: 0.2.8).

Main changes (see NEWS.md):

* Optional conditional models for one target variable: CART or penalised GAM
  for a continuous response (with training-only cross-validated engine
  selection), and a multinomial model for an unordered categorical response.
  These use the suggested packages 'rpart', 'mgcv' and 'nnet' only when
  requested, and the base package 'splines'.
* Unordered categorical columns are now placed on their latent copula
  coordinate in a data-driven order, so results for such columns differ from
  earlier versions; results for other data are unchanged.
* An algebraically equivalent, faster evaluation of the Student-t copula
  degrees-of-freedom likelihood.

## Test environments

* Local: Ubuntu Linux, R 4.3.3, `R CMD check --as-cran`.

## R CMD check results

0 errors | 0 warnings | 1 note

* Packages suggested but not available for checking: 'arrow', 'fst',
  'nanoparquet'. This NOTE is specific to the local check machine, where these
  optional I/O backends are not installed; the corresponding tests are skipped
  when a backend is unavailable.

## Reverse dependencies

There are currently no reverse dependencies.
