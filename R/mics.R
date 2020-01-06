#' Cross-spectrum based on mutual information
#' @export

mics <- function(
  spectro1= NULL, spectro2= NULL,
  feature1=  NULL, feature2= NULL,
  plot= F, zerofill= T, lagw = NULL,
  ...
){
  if (!is.null(spectro1) & !is.null(spectro2)){
    x <- spectro1
    y <- spectro2

    # Select longer signal as the static one
    if (ncol(x) >= ncol(y)) {
      static <- x
      lagged <- y
    }else {
      static <- y
      lagged <- x
    }
    # Generate band of indices for the lag series
    lagged <- cbind(lagged, rep(0, nrow(lagged)))
    band <- c(
      rep(ncol(lagged), ncol(static) - 1),
      1:(ncol(lagged) - 1),
      rep(ncol(lagged), ncol(static) - 1)
    )
    imax <- length(band) - ncol(static) + 1
    lags <- (ncol(static) - 1):(ncol(static) - imax)

    if (!is.null(lagw)){
      is <- (1:imax)[lags == lagw]
    }else{
      is <- 1:imax
    }

    lagstats <- lapply(is, function(i){
      iband <- (1:length(band) - 2 + i) %% length(band) + 1
      ycol <- band[iband][1:ncol(static)]
      if (zerofill){
        xcol <- 1:ncol(static)
        ycol <- band[iband][1:ncol(static)]
      }else{
        xcol <- (1:ncol(static))[-which(ycol == max(band))]
        ycol <- ycol[-which(ycol == max(band))]
      }
      x <- as.numeric(static[,xcol])
      y <- as.numeric(lagged[,ycol])
      mi_vector(x,y, ...)
    })
  }else if (!is.null(feature1) & !is.null(feature2)){
    x <- feature1
    y <- feature2

    # Select longer signal as the static one
    if (length(x) >= length(y)) {
      static <- x
      lagged <- y
    }else {
      static <- y
      lagged <- x
    }

    # Generate band of indices for the lag series
    lagged <- c(lagged, 0)
    band <- c(
      rep(length(lagged), length(static) - 1),
      1:(length(lagged) - 1),
      rep(length(lagged), length(static) - 1)
    )
    imax <- length(band) - length(static) + 1
    lags <- (length(static) - 1):(length(static) - imax)

    if (!is.null(lagw)){
      is <- (1:imax)[lags == lagw]
    }else{
      is <- 1:imax
    }

    lagstats <- lapply(is, function(i){
      iband <- (1:length(band) - 2 + i) %% length(band) + 1
      yis <- band[iband][1:length(static)]
      if (zerofill){
        xis <- 1:length(static)
        yis <- band[iband][1:length(static)]
      }else{
        xis <- (1:length(static))[-which(yis == max(band))]
        yis <- yis[-which(yis == max(band))]
      }
      x <- as.numeric(static[xis])
      y <- as.numeric(lagged[yis])
      mi_vector(x,y, ...)
    })
  }

  lagstats <- do.call("rbind", lagstats)
  if (!is.null(lagw)) {
    lagdf <- lagw
  }else{
    lagdf <- lags
  }
  lagstats <- data.frame(
    lag= lagdf,
    lagstats
  )

  if (plot){
    plot(mi ~ lag, data= lagstats, type= "l",
         ylab= "Mutual information", xlab= "lag (spectral windows)")
    plot(lc ~ lag, data= lagstats, type= "l",
         ylab= "Linear correlation", xlab= "lag (spectra windows")
  }


  return(lagstats)
}
