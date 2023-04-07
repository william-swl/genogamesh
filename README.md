
<!-- README.md is generated from README.Rmd. Please edit that file -->

# genogamesh

<!-- badges: start -->
<!-- badges: end -->

## installation

You can install the development version of `genogamesh` like so:

``` r
devtools::install_github("william-swl/genogamesh")
```

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
