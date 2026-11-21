# Get started with rthologs

rthologs looks up cross-species gene orthologs and returns a tidy R data
frame.

## Installation

Install the development version from GitHub with pak.

``` r

pak::pak("igordot/rthologs")
```

Load rthologs to run the examples.

``` r

library(rthologs)
```

## Orthologs

Use
[`orthologs()`](https://igordot.github.io/rthologs/reference/orthologs.md)
to find the orthologs of a set of genes. Give the genes, the source
species of the genes, and the target species.

``` r

odf <- orthologs("TP53", source = "human", target = "mouse")
odf
#>   input source_symbol target_symbol target_ncbi_id  target_ensembl_id
#> 1  TP53          TP53         Trp53          22059 ENSMUSG00000059552
```

The result is a standard data frame with one row for each known ortholog
of the input genes. The `target_` columns identify the ortholog.

The source and target can be any two species in the package.

``` r

odf <- orthologs("Trp53", source = "mouse", target = "rat")
odf
#>   input source_symbol target_symbol target_ncbi_id  target_ensembl_id
#> 1 Trp53         Trp53          Tp53          24842 ENSRNOG00000010756
```

### Input genes

The input can be gene symbols, NCBI IDs, or Ensembl IDs. Inputs are
case-sensitive and you can mix the input types in one vector.

``` r

genes <- c("TP53", "7157", "ENSG00000141510")
odf <- orthologs(genes, source = "human", target = "mouse")
odf
#>             input source_symbol target_symbol target_ncbi_id  target_ensembl_id
#> 1            TP53          TP53         Trp53          22059 ENSMUSG00000059552
#> 2            7157          TP53         Trp53          22059 ENSMUSG00000059552
#> 3 ENSG00000141510          TP53         Trp53          22059 ENSMUSG00000059552
```

An input with no direct match can match a gene synonym. To use only the
primary names, set `synonyms = FALSE`.

``` r

genes <- c("GAPDH", "CTGF", "SEPT7")
orthologs(genes, source = "human", target = "mouse")
#>   input source_symbol target_symbol target_ncbi_id  target_ensembl_id
#> 1 GAPDH         GAPDH         Gapdh          14433 ENSMUSG00000057666
#> 2  CTGF          CCN2          Ccn2          14219 ENSMUSG00000019997
#> 3 SEPT7       SEPTIN7       Septin7         235072 ENSMUSG00000001833
```

``` r

orthologs(genes, source = "human", target = "mouse", synonyms = FALSE)
#>   input source_symbol target_symbol target_ncbi_id  target_ensembl_id
#> 1 GAPDH         GAPDH         Gapdh          14433 ENSMUSG00000057666
```

### Output columns

Set `details = TRUE` to include the source gene IDs and the supporting
methods. The `support` and the `support_n`columns list the ortholog
sources and their count.

``` r

odf <- orthologs("TP53", source = "human", target = "mouse", details = TRUE)
dplyr::glimpse(odf)
#> Rows: 1
#> Columns: 11
#> $ input             <chr> "TP53"
#> $ source_taxon_id   <int> 9606
#> $ source_symbol     <chr> "TP53"
#> $ source_ncbi_id    <chr> "7157"
#> $ source_ensembl_id <chr> "ENSG00000141510"
#> $ target_taxon_id   <int> 10090
#> $ target_symbol     <chr> "Trp53"
#> $ target_ncbi_id    <chr> "22059"
#> $ target_ensembl_id <chr> "ENSMUSG00000059552"
#> $ support           <chr> "HGNC|Hieranoid|InParanoid|OMA|OrthoFinder|OrthoInsp…
#> $ support_n         <int> 9
```

### Multiple matches

An input can have more than one ortholog. Each ortholog gets its own
row.

``` r

orthologs(genes = c("Ptprc", "Hba-a1", "Gapdh"), source = "mouse", target = "zebrafish")
#>    input source_symbol target_symbol target_ncbi_id  target_ensembl_id
#> 1  Ptprc         Ptprc         ptprc         559154 ENSDARG00000071437
#> 2 Hba-a1        Hba-a1       hbae1.1          30597 ENSDARG00000089475
#> 3 Hba-a1        Hba-a1         hbaa2         497166 ENSDARG00000069735
#> 4  Gapdh         Gapdh         gapdh         317743 ENSDARG00000043457
```

By default, the result includes only the best matches. To get all the
potential orthologs, set `best = FALSE`.

``` r

odf <- orthologs(
  "Hba-a1",
  source = "mouse",
  target = "zebrafish",
  best = FALSE,
  details = TRUE
)
odf[, c("input", "target_symbol", "support_n")]
#>    input target_symbol support_n
#> 1 Hba-a1       hbae1.1         8
#> 2 Hba-a1         hbaa2         8
#> 3 Hba-a1         hbae4         4
```

Set `min_support` to keep pairs with a certain number of supporting
methods.

``` r

odf <- orthologs(
  "Hba-a1",
  source = "mouse",
  target = "zebrafish",
  best = FALSE,
  min_support = 6,
  details = TRUE
)
odf[, c("input", "target_symbol", "support_n")]
#>    input target_symbol support_n
#> 1 Hba-a1       hbae1.1         8
#> 2 Hba-a1         hbaa2         8
```

### Inputs without orthologs

By default, the output only includes genes with known orthologs. To keep
all inputs, set `orthologs_only = FALSE`.

``` r

orthologs(
  genes = c("Ptprc", "Cd34", "Fakegene"),
  source = "mouse",
  target = "zebrafish",
  orthologs_only = FALSE
)
#>      input source_symbol target_symbol target_ncbi_id  target_ensembl_id
#> 1    Ptprc         Ptprc         ptprc         559154 ENSDARG00000071437
#> 2     Cd34          <NA>          <NA>           <NA>               <NA>
#> 3 Fakegene          <NA>          <NA>           <NA>               <NA>
```

Genes without orthologs get one row with `NA` in all columns. `Cd34` is
a valid mouse gene, but the data contain no zebrafish ortholog for it.
`Fakegene` is not a mouse gene.

## Species

Use
[`species()`](https://igordot.github.io/rthologs/reference/species.md)
to list the available species and their accepted names.

``` r

species()
#>    taxon_id                           name          name_class
#> 1      6239         Caenorhabditis elegans     scientific name
#> 2      6239         Caenorhabditis elegans        species name
#> 3      7227        Drosophila melanogaster     scientific name
#> 4      7227        Drosophila melanogaster        species name
#> 5      7227                      fruit fly genbank common name
#> 6      7955                    Danio rerio     scientific name
#> 7      7955                    Danio rerio        species name
#> 8      7955                  leopard danio         common name
#> 9      7955                    zebra danio         common name
#> 10     7955                     zebra fish         common name
#> 11     7955                      zebrafish genbank common name
#> 12     8364             Xenopus tropicalis     scientific name
#> 13     8364             Xenopus tropicalis        species name
#> 14     8364           tropical clawed frog genbank common name
#> 15     8364            western clawed frog         common name
#> 16     9606                   Homo sapiens     scientific name
#> 17     9606                   Homo sapiens        species name
#> 18     9606                          human genbank common name
#> 19    10090                   Mus musculus     scientific name
#> 20    10090                   Mus musculus        species name
#> 21    10090                    house mouse genbank common name
#> 22    10090                          mouse         common name
#> 23    10116                     Norway rat genbank common name
#> 24    10116              Rattus norvegicus     scientific name
#> 25    10116              Rattus norvegicus        species name
#> 26    10116                      brown rat         common name
#> 27    10116                            rat         common name
#> 28    10116                           rats         common name
#> 29   559292       Saccharomyces cerevisiae        species name
#> 30   559292 Saccharomyces cerevisiae S288C     scientific name
```

## Data

The orthologs come from the [Alliance of Genome
Resources](https://www.alliancegenome.org/) (AGR). AGR combines
predictions from these sources: PANTHER, SonicParanoid, OrthoFinder,
Ensembl Compara, OrthoInspector, InParanoid, PhylomeDB, OMA, Hieranoid,
HGNC, and ZFIN.

The package uses the AGR stringent set that includes:

- Pairs that three or more methods support, with a best match in at
  least one direction.
- Pairs that two methods support, with a best match in both directions.
- All pairs from HGNC or ZFIN.

If you use the data, please cite: Alliance of Genome Resources
Consortium. Updates to the Alliance of Genome Resources central
infrastructure. *Genetics*. 2024. doi:
[10.1093/genetics/iyae049](https://doi.org/10.1093/genetics/iyae049).
