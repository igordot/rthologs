# Find gene orthologs

Find cross-species orthologs for a set of genes.

## Usage

``` r
orthologs(
  genes,
  source,
  target,
  synonyms = TRUE,
  best = TRUE,
  min_support = 1,
  orthologs_only = TRUE,
  details = FALSE
)
```

## Arguments

- genes:

  A character vector of gene symbols or NCBI/Ensembl gene IDs
  (case-sensitive).

- source:

  Species of `genes`, such as `"mouse"`. See
  [`species()`](https://igordot.github.io/rthologs/reference/species.md).

- target:

  Species in which to find orthologs.

- synonyms:

  If `TRUE` (default), check synonyms for input genes.

- best:

  If `TRUE` (default), keep only the best matches.

- min_support:

  Minimum number of sources that support each gene pair.

- orthologs_only:

  If `TRUE` (default), keep only matching pairs.

- details:

  If `TRUE`, return the source gene IDs and the support sources.

## Value

A data frame with one row per distinct input and matching gene pair.

## Details

The function uses each distinct input once. Different inputs for the
same gene each get their own rows. A valid gene can have no known
orthologs in the data.

The source AGR dataset includes orthologs from multiple sources, such as
PANTHER, SonicParanoid, and OrthoFinder. For pairs supported by three or
more methods, a best match in either direction is sufficient for the AGR
stringent set. For pairs supported by two methods, a best match in both
directions is required. The set also includes all pairs from ZFIN or
HGNC.

## References

Alliance of Genome Resources Consortium. Updates to the Alliance of
Genome Resources central infrastructure. Genetics. 2024. doi:
[doi:10.1093/genetics/iyae049](https://doi.org/10.1093/genetics/iyae049)
.

## Examples

``` r
orthologs(genes = "TP53", source = "human", target = "mouse")
#>   input source_symbol target_symbol target_ncbi_id  target_ensembl_id
#> 1  TP53          TP53         Trp53          22059 ENSMUSG00000059552
orthologs(genes = "Trp53", source = "mouse", target = "rat")
#>   input source_symbol target_symbol target_ncbi_id  target_ensembl_id
#> 1 Trp53         Trp53          Tp53          24842 ENSRNOG00000010756
orthologs(genes = c("Ptprc", "Cd34"), source = "mouse", target = "human")
#>   input source_symbol target_symbol target_ncbi_id target_ensembl_id
#> 1 Ptprc         Ptprc         PTPRC           5788   ENSG00000081237
#> 2  Cd34          Cd34          CD34            947   ENSG00000174059
```
