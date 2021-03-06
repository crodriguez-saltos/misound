% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chopAudio.R
\name{chopAudio}
\alias{chopAudio}
\title{Chop audio file into several segments}
\usage{
chopAudio(wave = NULL, file = NULL, prefix = NULL, timethrs,
  mindur = 0.01, plot = F, saveWav = T, outdir,
  buffersilence = NULL, silenceparams = NULL, ...)
}
\arguments{
\item{wave}{Wave file object}

\item{file}{Wave file address}

\item{prefix}{Prefix for audio output. This argument is required if file is
not given.}

\item{timethrs}{The time constant used to segment the audio.}

\item{mindur}{Minimum duration for segment to be considered valid.}

\item{outdir}{Directory where to save extracted segments.}

\item{buffersilence}{Name of audio file containing silence that is appended
at the beginning and end of each extracted segment of audio.}

\item{silenceparams}{Parameters for generating silence that will be appended
at the beginning and end of each extracted segment of audio.}

\item{...}{Arguments passed to seewave::timer.}
}
\description{
Chop audio file into several segments
}
\details{
Chops an audio file into several segments defined by a time
constant.

The script uses the seewave::timer function to ensure that the segment
boundaries do not cut across sounds.
}
\keyword{audio,segmentation}
