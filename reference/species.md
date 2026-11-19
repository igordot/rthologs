# List available species

List all the accepted names for the species in the package data.

## Usage

``` r
species(x = NULL)
```

## Arguments

- x:

  A species name or NCBI taxon ID. Use `NULL` (default) to list all
  species.

## Value

A data frame with `taxon_id`, `name`, and `name_class`.

## References

Schoch CL, Ciufo S, Domrachev M, et al. NCBI Taxonomy: a comprehensive
update on curation, resources and tools. Database (Oxford). 2020. doi:
[doi:10.1093/database/baaa062](https://doi.org/10.1093/database/baaa062)
.

## Examples

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
species("Mus musculus")
#>   taxon_id         name          name_class
#> 1    10090 Mus musculus     scientific name
#> 2    10090 Mus musculus        species name
#> 3    10090  house mouse genbank common name
#> 4    10090        mouse         common name
species("rat")
#>   taxon_id              name          name_class
#> 1    10116        Norway rat genbank common name
#> 2    10116 Rattus norvegicus     scientific name
#> 3    10116 Rattus norvegicus        species name
#> 4    10116         brown rat         common name
#> 5    10116               rat         common name
#> 6    10116              rats         common name
species(10090)
#>   taxon_id         name          name_class
#> 1    10090 Mus musculus     scientific name
#> 2    10090 Mus musculus        species name
#> 3    10090  house mouse genbank common name
#> 4    10090        mouse         common name
```
