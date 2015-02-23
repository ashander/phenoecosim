#' G matrix
#' @param Gaa additive genetic variance in RN intercept
#' @param Gbb additive genetic variance in RN slope
#' @param Gab additive genetic covariance between RN slope and height
#' @export
G <- function(Gaa, Gbb, Gab)
    matrix( c(Gaa, Gab, Gab, Gbb), nrow=2) #eqn 3a
