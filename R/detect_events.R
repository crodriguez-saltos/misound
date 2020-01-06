#'@export

detect_events <- function(wave= NULL, file= NULL,
                         mindur= 0.01, ...){
  # Detect sound events----
  if (is.null(wave)){
    sound <- tuneR::readWave(file)
  }else{
    sound <- wave
  }

  # Add silence at beginning and end----
  # This helps with detecting signals
  sound_dur <- seewave::duration(sound)
  silence_dur <- 0.2
  sound <- seewave::addsilw(wave = sound, at = "start",
                            d = silence_dur, output = "Wave")
  sound <- seewave::addsilw(wave = sound, at = "end",
                            d = silence_dur, output = "Wave")

  segments <- seewave::timer(wave = sound, plot= F, dmin= mindur, ...)

  # Ensure that starts always precede ends----
  starts.first <- segments$s.start[1]
  starts.last <- segments$s.start[length(segments$s.start)]
  ends.first <- segments$s.end[1]
  ends.last <- segments$s.end[length(segments$s.end)]

  if (starts.first > ends.first){
    segments$s.start <- c(0, segments$s.start)
  }

  if (starts.last > ends.last){
    segments$s.end <- c(segments$s.end, duration(sound))
  }

  # Extract timestamps----
  timestamps <- data.frame(start= segments$s.start,
                           end= segments$s.start + segments$s)
  timestamps <- timestamps - silence_dur

  timestamps[timestamps < 0] <- 0
  timestamps[timestamps > sound_dur] <- sound_dur

  return(timestamps)
}
