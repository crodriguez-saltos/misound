# misound

`misound` is a R package that allows comparing sounds using mutual information. It includes functions for batch analysis of syllables.

## Installing misound
`misound` can be installed directly from GitHub. The packages `devtools`, `knitr`, `rmarkdown`, and `seewave` are required. If any of those packages is not present in R, install it with the following code:

```{r}
install.packages("devtools")
install.packages("knitr")
install.packages("rmarkdown")
install.packages("seewave")
```

Then install `misound` with the following code:

```{r}
install_github(repo= "crodriguez-saltos/misound", build_vignettes= TRUE)
```

## Getting started
The package comes with a useful vignette for getting started. Check the vignette by using the following code:

```{r}
browseVignettes(package= "misound")
```
