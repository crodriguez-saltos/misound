#' Short-term FFT only for signal portions of a sound

#' @export

signal_spectro <- function(
  sound,
  plot= F,
  segment_params = list(wl= 512, ovlp= 0, threshold= 5),
  ...
){
  # Get spectrogram
  spectro <- seewave::spectro(
    wave = sound,
    dB = NULL,
    plot= F,
    wl= segment_params$wl,
    ovlp= segment_params$ovlp,
    ...)

  # Block silences
  if (segment_params$threshold != 0){
    signal <- detect_events(
      wave= sound,
      msmooth= c(segment_params$wl, segment_params$ovlp)
    )
  }else{
    signal <- data.frame(start= 0, end= duration(sound))
  }

  whichstart <- cut(signal$start, breaks = spectro$time)
  if (signal$start[1] == 0){
    whichstart[1] <- levels(whichstart)[1]
  }
  if (round(seewave::duration(sound), 6) == round(tail(signal$end, n= 1), 6)){
    signal$end[length(signal$end)] <- tail(spectro$time, n= 1)
  }
  whichend <- cut(signal$end, breaks = spectro$time)
  signal <- rbind(as.numeric(whichstart), as.numeric(whichend))
  signal <- apply(signal, 2, function(x) seq(x[1], x[2]))
  if (class(signal) == "matrix"){
    signal <- as.data.frame(signal)
    signal <- as.list(signal)
  }
  signal <- do.call("c", signal)

  spectro.signal <- lapply(1:length(spectro$time), function(i){
      if (is.element(i, signal)){
        spectro$amp[,i]
      }else{
        rep(NA, length(spectro$freq))
      }
  })
  spectro.signal <- as.matrix(do.call("cbind", spectro.signal))

  spectro.signal <- list(
    time= spectro$time,
    freq= spectro$freq,
    amp= spectro.signal
  )

  if (plot){
    image(t(spectro.signal$amp))
  }
  return(spectro.signal)

}
