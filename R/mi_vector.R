#' Calculate mutual information between vectors
#'
#' @param x,y Vectors to be compared using mutual information.
#' @param mitype Type of mutual information estimator to be used. See details.
#' @param minorm Should mutual information be normalized? See details.
#' @param autoswitch If KDE fails, should estimator be switched to Jackknife?.
#'   See details.
#' @param lc Should linear correlation between x and y be calculated?
#'
#' @details Two types of estimators of mutual information can be chosen by the
#'   user: KDE and Jackknife. The latter tends to give consistent results, but
#'   it is slower.
#'
#'   When normalization is applied, mutual information is divided by the average
#'   of the entropies of signals x and y.
#'
#'   Sometimes KDE fails to produce a mutual information estimate. A common case
#'   is when having a sample size too low for KDE estimation.
#'
#'   Results of mutual information are expressed in nats.
#'
#' @export

mi_vector <- function(x,y, mitype= "kde",
                      minorm= F, autoswitch= F, lc= T){
  if (lc){
    cor <- cor(x, y)
  }else{
    cor <- NA
  }

  if (mitype == "kde"){
    mi <- tryCatch(expr = {
      mikde(x, y)
    }, error= function(e){
      return(NA)
    })
    if (is.infinite(mi) & autoswitch){
      mi <- JMI::JMI(x, y, BN = 10)$mi
      mitype <- "jacknife"
      print("Not suitable for KDE. Switching to Jackknife")
    }
  }else if (mitype == "jackknife"){
    mi <- JMI::JMI(x, y, BN = 10)$mi
  }else if (mitype == "none"){
    mi <- NA
    minorm <- F
  }
  if (minorm){
    if (mitype == "kde"){
      h1 <- tryCatch(expr = {
        mikde(x, x)
      }, error= function(e){
        return(NA)
      })
      h2 <- tryCatch(expr = {
        mikde(y, y)
      }, error= function(e){
        return(NA)
      })
    }else if (mitype == "jackknife"){
      h1 <- JMI::JMI(x, x, BN = 10)$mi
      h2 <- JMI::JMI(y, y, BN = 10)$mi
    }
    minorm <- 2 * mi / (h1 + h2)
  }else{
    minorm <- NA
  }

  nlc <- ifelse(
    test= is.infinite(mi) | is.na(mi),
    yes = NA,
    no = sqrt(1 - exp(-2 * mi))
  )
  return(data.frame(mi= mi, nlc= nlc, lc= cor, minorm= minorm, mitype= mitype))
}
