# rthologs: Gene Orthologs for Model Organisms

rthologs maps genes to their orthologs in another species.

- Supports common model organisms, including human, mouse, zebrafish,
  and fly.
- Looks up genes by symbol, NCBI (Entrez) ID, or Ensembl ID.
- Returns matches as a tidy R data frame you can filter, summarize, and
  join to your own data.
- Works offline.

## Installation

Install the development version from GitHub.

``` r

pak::pak("igordot/rthologs")
```

## Usage

Find the mouse ortholog of the human gene `TP53`.

``` r

rthologs::orthologs("TP53", source = "human", target = "mouse")
```

Use
[`rthologs::species()`](https://igordot.github.io/rthologs/reference/species.md)
to list the accepted species names.

See the [package documentation](https://igordot.github.io/rthologs/) for
more details.
