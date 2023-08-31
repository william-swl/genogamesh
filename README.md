
<!-- README.md is generated from README.Rmd. Please edit that file -->

# genogamesh

<!-- badges: start -->
<!-- badges: end -->

## installation

You can install the development version of `genogamesh` like so:

``` r
devtools::install_github("william-swl/genogamesh")
```

## parse bioinfomatic data

- parse the output of SingleR

``` r
# SingleR(test, ref) %>% parse_SingleR()
```

- parse somatic hypermutation from igblast output

``` r
# parse_IgBlast_shm('igblast_out.txt')
```

- parse sequences from CellRanger vdj output

``` r
# parse_CellRanger_vdjseq(df)
# parse_CellRanger_vdjseq(df, file='seq.csv')
# parse_CellRanger_vdjseq(df, file='seq.fa', fa_content='seq_orf_nt')
```

- parse sequences from ANARCI vdj output

``` r
# parse_ANARCI_aaseq(df, chain='H')
# parse_ANARCI_aaseq(df, chain='L')

# keep the ab numbering
# parse_ANARCI_aaseq(df, chain='H', keep_number=TRUE)
```

- parse vcf file with the help of reference genome and annotations. It
  is still under development, which can not process more than 3 nt
  substitutions in a single record row of vcf file, and can not process
  indels

``` r
# vcf <- read_vcf(...)
# fa <- read_fa(..)
# gff <- read_gff(...)
# parse_vcf(vcf, fa, gff)
```

## shortcuts for bioinfomatic pipelines

- add SingleR celltype annotation for Seurat object

``` r
# SE <- SingleR_SE(SE, SEref)
```

- reduction from raw Seurat object created by read count matrix,
  including normalization, variable features calling, scaling, PCA and
  UMAP

``` r
# SE <- reduction_SE(SE)
```

- translate nucleotides into amino acids from the first character

``` r
nt2aa(c("ATGAAA", "TTGCCC", "CTGTTT"))
#> [1] "MK" "LP" "LF"
```

- build antigen map from sera titer data

``` r
# antigen_map(data, sera_meta, ag_meta, seed=14)
```

## IO

    read_fasta()

    write_fasta()

    read_vcf()

    read_gff()

## S4 classes in `genogamesh`

### mutstr

- a S4 class to manipulate mutation strings
- support set operations

``` r
raw_mut_string <- c(
  variant1 = "T10I,D20N,Q30E,A40T,P50L,G60R",
  variant2 = "T10I,D20-,Q30E,A40T,P50L,G60R,S80R",
  variant3 = "T10A,D20G,Q30E,A40T,P50L,G60R"
)

m <- mutstr(raw_mut_string, sep = ",")

m
#> mutstr 3
#>   @ names: variant1 variant2 variant3
#>   @ sep: ,
#>   @ mstr:
#>     [1] T10I,D20N,Q30E,A40T,P50L,G60R
#>     [2] T10I,D20-,Q30E,A40T,P50L,G60R,S80R
#>     [3] T10A,D20G,Q30E,A40T,P50L,G60R
#>   @ mut:
#>     [1] T10I D20N Q30E A40T P50L G60R
#>     [2] T10I D20- Q30E A40T P50L G60R S80R
#>     [3] T10A D20G Q30E A40T P50L G60R

names(m)
#> [1] "variant1" "variant2" "variant3"

mstr(m)
#>                             variant1                             variant2 
#>      "T10I,D20N,Q30E,A40T,P50L,G60R" "T10I,D20-,Q30E,A40T,P50L,G60R,S80R" 
#>                             variant3 
#>      "T10A,D20G,Q30E,A40T,P50L,G60R"

mut(m)
#> $variant1
#> [1] "T10I" "D20N" "Q30E" "A40T" "P50L" "G60R"
#> 
#> $variant2
#> [1] "T10I" "D20-" "Q30E" "A40T" "P50L" "G60R" "S80R"
#> 
#> $variant3
#> [1] "T10A" "D20G" "Q30E" "A40T" "P50L" "G60R"

m[1:2]
#> mutstr 2
#>   @ names: variant1 variant2
#>   @ sep: ,
#>   @ mstr:
#>     [1] T10I,D20N,Q30E,A40T,P50L,G60R
#>     [2] T10I,D20-,Q30E,A40T,P50L,G60R,S80R
#>   @ mut:
#>     [1] T10I D20N Q30E A40T P50L G60R
#>     [2] T10I D20- Q30E A40T P50L G60R S80R

m[[2]]
#> [1] "T10I" "D20-" "Q30E" "A40T" "P50L" "G60R" "S80R"

intersect(m, m[1])
#> mutstr 3
#>   @ names: variant1 variant2 variant3
#>   @ sep: ,
#>   @ mstr:
#>     [1] T10I,D20N,Q30E,A40T,P50L,G60R
#>     [2] T10I,Q30E,A40T,P50L,G60R
#>     [3] Q30E,A40T,P50L,G60R
#>   @ mut:
#>     [1] T10I D20N Q30E A40T P50L G60R
#>     [2] T10I Q30E A40T P50L G60R
#>     [3] Q30E A40T P50L G60R

setdiff(m, m[1])
#> mutstr 3
#>   @ names: variant1 variant2 variant3
#>   @ sep: ,
#>   @ mstr:
#>     [1] 
#>     [2] D20-,S80R
#>     [3] T10A,D20G
#>   @ mut:
#>     [1] 
#>     [2] D20- S80R
#>     [3] T10A D20G

union(m, m[1])
#> mutstr 3
#>   @ names: variant1 variant2 variant3
#>   @ sep: ,
#>   @ mstr:
#>     [1] T10I,D20N,Q30E,A40T,P50L,G60R
#>     [2] T10I,D20-,Q30E,A40T,P50L,G60R,S80R,D20N
#>     [3] T10A,D20G,Q30E,A40T,P50L,G60R,T10I,D20N
#>   @ mut:
#>     [1] T10I D20N Q30E A40T P50L G60R
#>     [2] T10I D20- Q30E A40T P50L G60R S80R D20N
#>     [3] T10A D20G Q30E A40T P50L G60R T10I D20N
```
