#' Take similarity measurements in batch
#' @export

similarity_batch <- function(
  sound1, sound2, s1.sylls, s2.sylls,
  spectroargs= list(wl= 512, ovlp= 0, wn= "hanning"),
  alignment= "lc", datatype= "spectro"
    ){

  # spectrogram function
  spectrobatch <- function(x,...){
    seewave::spectro(dB= NULL, plot= F, ...)
  }

  featurebatch <- function(x...){
    spectral_features(
      ...,
      windowargs = list(
        wl= spectroargs$wl, ovlp= spectroargs$ovlp, threshold= 0
      ))
  }

  sylls.comp <- expand.grid(1:nrow(s1.sylls), 1:nrow(s2.sylls))
  simstats <- lapply(1:nrow(sylls.comp), function(i){
    tryCatch(
      expr = {
        print(
          paste("Scoring similarity for comparison", i, "of", nrow(sylls.comp))
        )
        x <- as.numeric(sylls.comp[i,])

        wave1 <- cutw(
          sound1,
          from= s1.sylls$start[x[1]],
          to= s1.sylls$end[x[1]],
          output= "Wave"
        )

        wave2 <- cutw(
          sound2,
          from= s2.sylls$start[x[2]],
          to= s2.sylls$end[x[2]],
          output= "Wave"
        )

        if (datatype == "spectro"){
          s1syll <- do.call(
            what = "spectrobatch",
            args = c(
              spectroargs,
              wave= wave1
            )
          )$amp

          s2syll <- do.call(
            what = "spectrobatch",
            args = c(
              spectroargs,
              wave= wave2
            )
          )$amp

          # Normalize spectra
          s1syll <- (s1syll - min(s1syll)) / (max(s1syll) - min(s1syll))
          s2syll <- (s2syll - min(s2syll)) / (max(s2syll) - min(s2syll))

          # Mutual information
          if (alignment == "lc"){
            sims <- mics(spectro1 = s1syll, spectro2 = s2syll, mitype = "none")
            lclag <- T
            milag <- F
          }else if (alignment == "mi"){
            sims <- mics(spectro1 = s1syll, spectro2 = s2syll, mitype = "kde")
            lclag <- F
            milag <- T
          }else if (alignment == "both"){
            sims <- mics(spectro1 = s1syll, spectro2 = s2syll, mitype = "kde")
            lclag <- T
            milag <- T
          }
        }else if (datatype == "h"){
          feat1 <- seewave::csh(
            wave = wave1,
            wl = spectroargs$wl,
            ovlp = spectroargs$ovlp,
            plot= F)[,2]

          feat2 <- seewave::csh(
            wave = wave2,
            wl = spectroargs$wl,
            ovlp = spectroargs$ovlp,
            plot= F)[,2]

          # Mutual information
          if (alignment == "lc"){
            sims <- mics(feature1 = feat1, feature2= feat2, mitype = "none")
          }else if (alignment == "mi"){
            sims <- mics(feature1 = feat1, feature2= feat2, mitype = "kde",
                         minorm= T)
          }else if (alignment == "both"){
            sims <- mics(feature1 = feat1, feature2= feat2, mitype = "kde",
                         minorm= T)
          }
        }

        if (alignment == "lc"){
          lclag <- T
          milag <- F
        }else if (alignment == "mi"){
          lclag <- F
          milag <- T
        }else if (alignment == "both"){
          lclag <- T
          milag <- T
        }

        if (lclag){
          maxlc <- max(sims$lc[is.finite(sims$lc)], na.rm= T)
          maxlcdf <- sims[which(sims$lc == maxlc),]
          maxlcdf <- data.frame(comparison= i, maxlcdf, maxtype= "lc")
        }

        if (milag){
          maxmi <- max(sims$mi[is.finite(sims$mi)], na.rm= T)
          maxmidf <- sims[which(sims$mi == maxmi),]
          maxmidf <- data.frame(comparison= i, maxmidf, maxtype= "mi")
        }

        if (alignment == "lc"){
          maxdf <- maxlcdf
        }else if (alignment == "mi"){
          maxdf <- maxmidf
        }else if (alignment == "both"){
          maxdf <- rbind(maxlcdf, maxmidf)
        }

        if (alignment == "lc"){
          if (datatype == "spectro"){
            maxdf <- mics(spectro= s1syll, spectro2= s2syll, mitype= "kde",
                          lag= maxdf$lag, minorm= T)
            maxdf$maxwindows <- max(ncol(s1syll), ncol(s2syll))
          }else if (datatype == "h"){
            maxdf <- mics(feature1 = feat1, feature2= feat2, mitype= "kde",
                          lag= maxdf$lag, minorm= T)
            maxdf$maxwindows <- max(length(feat1), length(feat2))
          }
          maxdf <- data.frame(comparison= i, maxdf, maxtype= "lc")
        }

        maxdf$s1syll <- sylls.comp[i,1]
        maxdf$s2syll <- sylls.comp[i,2]

        return(maxdf)
      }, error= function(e){
        NA
      }
    )
  })
  simstats <- do.call("rbind", simstats)
  return(simstats)
}
