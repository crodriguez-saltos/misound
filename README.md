# misound

`misound` is a R package that allows comparing sounds using mutual information. It includes functions for batch analysis of syllables.

## Installing misound
`misound` can be installed directly from GitHub. The package `devtools` is required. If that package is not present in R, install it with the following line of code:

```{r}
install.packages("devtools")
```

Then install `misound` with the following code:

```{r}
install_github(repo= "crodriguez-saltos/misound", build_vignettes= TRUE)
```

## Getting started
The package comes a useful vignettes for getting started. Check the vignettes by using the following code:

```{r}
browseVignettes(package= "misound")
```
