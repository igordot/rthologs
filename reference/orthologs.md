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
  orthologs_only = TRUE
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

  If `TRUE` (default), check gene synonyms for inputs with no direct
  match.

- best:

  If `TRUE` (default), keep only the best matches, which can differ
  between search directions.

- min_support:

  Minimum number of methods or sources that support a gene pair. Use one
  positive whole number. The maximum possible count depends on the
  species pair.

- orthologs_only:

  If `TRUE` (default), keep only inputs with at least one matching pair.

## Value

A data frame with one row per distinct input and matching gene pair, in
input order. With `orthologs_only = FALSE`, each distinct nonmissing
input with no matching pair gets one row with its value in `input` and
`NA` in all other columns.

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
#>   input source_taxon_id source_symbol source_ncbi_id source_ensembl_id
#> 1  TP53            9606          TP53           7157   ENSG00000141510
#>   target_taxon_id target_symbol target_ncbi_id  target_ensembl_id
#> 1           10090         Trp53          22059 ENSMUSG00000059552
#>                                                                                    support
#> 1 HGNC|Hieranoid|InParanoid|OMA|OrthoFinder|OrthoInspector|PANTHER|PhylomeDB|SonicParanoid
#>   support_n
#> 1         9
orthologs(genes = c("Ptprc", "Cd34"), source = "mouse", target = "human")
#>   input source_taxon_id source_symbol source_ncbi_id  source_ensembl_id
#> 1 Ptprc           10090         Ptprc          19264 ENSMUSG00000026395
#> 2  Cd34           10090          Cd34          12490 ENSMUSG00000016494
#>   target_taxon_id target_symbol target_ncbi_id target_ensembl_id
#> 1            9606         PTPRC           5788   ENSG00000081237
#> 2            9606          CD34            947   ENSG00000174059
#>                                                                                                    support
#> 1 Ensembl Compara|HGNC|Hieranoid|InParanoid|OMA|OrthoFinder|OrthoInspector|PANTHER|PhylomeDB|SonicParanoid
#> 2 Ensembl Compara|HGNC|Hieranoid|InParanoid|OMA|OrthoFinder|OrthoInspector|PANTHER|PhylomeDB|SonicParanoid
#>   support_n
#> 1        10
#> 2        10
orthologs(genes = "Trp53", source = "mouse", target = "rat")
#>   input source_taxon_id source_symbol source_ncbi_id  source_ensembl_id
#> 1 Trp53           10090         Trp53          22059 ENSMUSG00000059552
#>   target_taxon_id target_symbol target_ncbi_id  target_ensembl_id
#> 1           10116          Tp53          24842 ENSRNOG00000010756
#>                                                                                               support
#> 1 Ensembl Compara|Hieranoid|InParanoid|OMA|OrthoFinder|OrthoInspector|PANTHER|PhylomeDB|SonicParanoid
#>   support_n
#> 1         9
```
