# rthologs: Gene Orthologs for Model Organisms

<!-- badges: start -->
[![R-CMD-check](https://github.com/igordot/rthologs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/igordot/rthologs/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/igordot/rthologs/graph/badge.svg)](https://app.codecov.io/gh/igordot/rthologs)
<!-- badges: end -->

rthologs maps genes to their orthologs in other species.

* Supports common model organisms, including human, mouse, zebrafish, and fly.
* Looks up genes by symbol, NCBI (Entrez) ID, or Ensembl ID.
* Returns matches as a tidy R data frame you can filter, summarize, and join to your own data.
* Works offline.

## Installation

Install the development version from GitHub.

```r
pak::pak("igordot/rthologs")
```

## Usage

Find the mouse ortholog of the human gene `TP53`.

```r
rthologs::orthologs("TP53", source = "human", target = "mouse")
```

Use `rthologs::species()` to list the accepted species names.

See the [vignette](https://igordot.github.io/rthologs/articles/rthologs.html) for more examples.
