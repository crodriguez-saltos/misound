#' Compare two sounds based using mutual information
#'
#' @param sound1 Sound 1. Wave object.
#' @param sound2 Sound 2. Wave object.
#' @param features Features data frame.
#' @param miN Number of neighbors for estimation of mutual information
#' @param d1 Discretized spectrogram of sound 1.
#' @param d2 Discretized spectrogram of sound 2.
#' @param norm Normalize mutual information?
#' @param align Should time traces be aligned?
#' @param type Type of mutual information estimator. See details.
#' @param nlcc Should the non-linear correlation coefficient be reported?. See details.
#' @param plot Plot
#' @param ... Additional arguments passed to infotheo::mutinformation.
#'
#' @details
#' The types of estimation of mutual information are:
#'
#' symba -  Estimate mutual information using seewave::symba.
#' infotheo - Estimate mutual information using package infotheo.
#'
#' The non-linear correlation coefficient is that of XXXX.
#'
#'
#'@export

sim_MI <- function(
  sound1= NULL,
  sound2= NULL,
  features= NULL,
  miN= 10,
  d1= NULL,
  d2= NULL,
  type= "jackknife",
  nlcc = F,
  plot= F,
  ...){

  # Import sound and get power spectra
  a <-tuneR::readWave(sound1)
  b <-tuneR::readWave(sound2)

  a.spec <- seewave::spectro(a, plot= F)$amp
  b.spec <- seewave::spectro(b, plot= F)$amp

  if (type == "symba"){
    mi <- matrix(data = NA, nrow = ncol(spectro1), ncol = ncol(spectro2))
    h1 <- mi
    h2 <- mi
    nmi <- rep(NaN, ncol(spectro1) * ncol(spectro2))
    for (i in 1:ncol(spectro1)){
      for (j in 1:ncol(spectro2)){
        pixel <- symba(x = spectro1[,i], y = spectro2[,j], plot = F)
        mi[i,j] <- pixel$I
        h1[i,j] <- pixel$h1
        h2[i,j] <- pixel$h2
        nmi[i,j] <- pixel$I / pixel$h1
      }
    }
  }else if (type == "infotheo"){
    # Mark the columns in the spectra that correspond to signal
    sound1.signal <- apply(sound1, 2, function(x) !all(is.na(x)))
    sound2.signal <- apply(sound2, 2, function(x) !all(is.na(x)))

    # Remove silence from signals
    sound1.filtered <- sound1[,sound1.signal]
    sound2.filtered <- sound2[,sound2.signal]

    # Estimate mutual information on signal
    mi <- infotheo::mutinformation(
      as.data.frame(cbind(sound1.filtered, sound2.filtered)),
      ...
    )

    # Select comparisons of interest
    sound1sound2 <- mi[
      1:ncol(sound1.filtered),
      (ncol(sound1.filtered) + 1):(ncol(sound1.filtered) + ncol(sound2.filtered))
      ]

    # Find the coordinates of the output matrix that correspond to each
    # syllable per syllable comparison
    is.onset <- function(x){
      c(x[1], !x[1:(length(x) - 1)] & x[2:length(x)])
    }
    is.offset <- function(x){
      c(x[1:(length(x) - 1)] & !x[2:length(x)], x[length(x)])
    }
    sound1.sylls <- data.frame(
      on = which(is.onset(sound1.signal)),
      off = which(is.offset(sound1.signal))
    )
    sound2.sylls <- data.frame(
      on = which(is.onset(sound2.signal)),
      off = which(is.offset(sound2.signal))
    )
    sylls.comp <- expand.grid(1:nrow(sound1.sylls), 1:nrow(sound2.sylls))
    sylls.coord <- mapply(FUN = function(x, y){
      x[y,]
    },
    x= list(sound1= sound1.sylls, sound2= sound2.sylls),
    y= sylls.comp
    )
    sylls.coord <- do.call("data.frame", sylls.coord)
    colnames(sylls.coord) <- c("sound1.on", "sound1.off", "sound2.on", "sound2.off")

    # Map out the corresponding coordinates in the mutual information matrix
    dursil <- function(x){
      x[,2] - x[,1]
    }
    sound1.dur <- dursil(sound1.sylls)
    sound2.dur <- dursil(sound2.sylls)
    sound1.sylls.mi <- data.frame(
      on= cumsum(c(1, sound1.dur[1:(length(sound1.dur) - 1)] + 1))
    )
    sound1.sylls.mi$off <- sound1.sylls.mi + sound1.dur
    sound2.sylls.mi <- data.frame(
      on= cumsum(c(1, sound2.dur[1:(length(sound2.dur) - 1)] + 1))
    )
    sound2.sylls.mi$off <- sound2.sylls.mi + sound2.dur

    sylls.comp.mi <- expand.grid(1:nrow(sound1.sylls.mi), 1:nrow(sound2.sylls.mi))
    sylls.coord.mi <- mapply(FUN = function(x, y){
      x[y,]
    },
    x= list(sound1= sound1.sylls.mi, sound2= sound2.sylls.mi),
    y= sylls.comp.mi
    )
    sylls.coord.mi <- do.call("data.frame", sylls.coord.mi)
    colnames(sylls.coord.mi) <- c("sound1.on", "sound1.off", "sound2.on", "sound2.off")

    # Fill the output matrix with MI values
    totalbins <- ncol(sound1) * ncol(sound2)
    mi <- rep(NA, totalbins)
    mi <- matrix(mi, nrow=ncol(sound1))

    for (i in 1:nrow(sylls.coord)){
      mi[
        sylls.coord$sound1.on[i]:sylls.coord$sound1.off[i],
        sylls.coord$sound2.on[i]:sylls.coord$sound2.off[i]
        ] <- sound1sound2[
          sylls.coord.mi$sound1.on[i]:sylls.coord.mi$sound1.off[i],
          sylls.coord.mi$sound2.on[i]:sylls.coord.mi$sound2.off[i]
          ]
    }
  }else if (type == "ksg"){
    if(plot){
      # Spectrograms
      bioacoustics::spectro(wave = a)
      bioacoustics::spectro(wave = b)
      plot(a.spec, b.spec)

      dd <- data.frame(
        sound1= as.numeric(a.spec),
        sound2= as.numeric(b.spec)
      )
      p <- ggplot(dd, aes(x=sound1, y=sound2) ) +
        geom_bin2d() +
        theme_bw()
      print(p)

      plot(a.spec, b.spec)
    }

    # Mutual information
    mi <- FNN::mutinfo(a.spec, b.spec, k = miN)
    #h1 <- FNN::mutinfo(a.spec, a.spec, k = miN)

    # Correlation
    #cor <- cor(as.numeric(a.spec), as.numeric(b.spec))

    # Output
    #print("Correlation coefficient is:")
    #print(cor)
    #print("Mutual information is")
    #print(mi)
    #print("Normalized mutual information is")
    #nmi <- mi/h1
    #print(mi/h1)
    #return(data.frame(corr= cor, nmi= nmi, mi= mi))
  }else if (type == "kde"){
    xy <- MASS::kde2d(x = x, y = y, n = 100)
    pxy <- xy$z
    pxy <- (pxy - min(pxy)) / (max(pxy) - min(pxy))
    pxy <- pxy / sum(pxy)

    px <- replicate(ncol(pxy), rowSums(pxy))
    py <- t(replicate(nrow(pxy), colSums(pxy)))

    if(plot){
      image(pxy)
    }

    mi <- sum(pxy * log(pxy / (px * py)), na.rm = T)
  }else if (type == "jackknife"){
    mi <- JMI::JMI(x, y)$mi
  }

  if (nlcc){
    nlcc <- sqrt(1 - exp(-2 * mi))
    return(data.frame(mi= mi, nlcc= nlcc))
  }else{
    return(mi)
  }
}
