#' Get spectral features of a wave file
#'
#' @param wave Wave object.
#' @param thresh Threshold for separating signal from noise. % of maximum
#'   amplitude.
#' @param ... Other \link[misound]{signal_spectro} arguments.
#'
#' @export

spectral_features <- function(
  wave,
  windowargs= list(wl= 512, ovlp= 0, threshold= 5),
  ...
  ){

  # Get spectrogram
  spectro <- signal_spectro(
    sound= wave,
    segment_params = windowargs,
    ...
  )

  # General features
  props <- apply(spectro$amp, 2, function(x){
    if (!all(is.na(x))){
      seewave::specprop(spec= x, f= wave@samp.rate)
    }else{
      list(status= NA)
    }
  })

  for (i in 1:length(props)){
    props[[i]]$time <- spectro$time[i]
  }

  # Fundamental frequency
  fundfun <- function(...){
    seewave::fund(
      wave, plot = F, ...
    )
  }

  funds <- do.call("fundfun", windowargs)
  for (i in 1:length(props)){
    if (!(names(props[[i]])[1] == "status")){
      props[[i]]$fund <- funds[i,2]
    }
  }

  # Spectral width
  specw <- function(spec, f) {
    sum(spec[, 2] * (spec[, 1] - f)^2)/sum(spec[, 2])
  }

  for (i in 1:length(props)){
    if (!(names(props[[i]])[1] == "status")){
      specmat <- cbind(spectro$freq, spectro$amp[,i])
      props[[i]]$fund_w <- specw(specmat, props[[i]]$fund[i])
      props[[i]]$mode_w <- specw(specmat, props[[i]]$mode/1000)
    }
  }

  # Consolidate measurements into one data frame
  props.eliminate <- sapply(props, function(x){
    names(x)[1] == "status"
  })
  props <- props[!props.eliminate]
  props <- lapply(props, as.data.frame)
  props <- do.call("rbind", props)
  eliminatel <- length(which(props.eliminate))
  nfullrows <- nrow(props)
  props[nrow(props) + eliminatel,] <- NA
  props$time[(nfullrows + 1):nrow(props)] <- spectro$time[props.eliminate]
  props <- props[order(props$time),]
  timecol <- which(colnames(props) == "time")
  props <- data.frame(time= props$time, props[,-timecol])
  rownames(props) <- 1:nrow(props)

  return(props)
}
