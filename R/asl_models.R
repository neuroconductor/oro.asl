#' @title Kinetic Models for Arterial Spin Labeling
#' 
#' @description Theoretical models for the acquisition of pulsed ASL data.
#' 
#' See Buxton etal. (1998) for more details.
#' 
#' @rdname asl2p
#' @aliases asl2p asl3p
#' @param beta is a vector (lenth = 2 or 3) of kinetic parameters.
#' @param TI is a vector of inversion times (TIs) in seconds.
#' @param T1b is the longitudinal relaxation time of blood.
#' @param T1 is the longitudinal relaxation time of tissue.
#' @param tau is the first TI in the second paper (?).
#' @param alpha is the fraction of maximum possible change in the longitudinal
#' magnetization that was achieved.
#' @param lambda is the equilibrium tissue/blood partition coefficient of
#' water.
#' @param Mob is the equilibrium magnetization of arterial blood.
#' @return A theoretical curve of pulsed ASL signal as a function of time.
#' @author Brandon Whitcher \email{bjw34032@@users.sourceforge.net}
#' @references Buxton, R.B., Frank, L.R., Wong, E.C., Siewert, B., Warach, S.
#' and Edelman, R.R. (1998) A General Kinetic Model for Quantitative Perfusion
#' Imaging with Arterial Spin Labeling, \emph{Magnetic Resonance in Medicine},
#' \bold{40}, 383-396.
#' @keywords models
#' @examples
#' TI <- seq(0, 3, length=256)
#' beta <- c(f=0.8, deltaT=0.5)
#' plot(TI, asl2p(beta, TI), type="l", lwd=2,
#'      xlab="Time (sec)", ylab="Signal", main="Pulsed ASL: Standard Model")
#' abline(v=1.5, lwd=2, col=2)
#' @export
asl2p <- function(beta, TI, T1b=1.3, T1=1.0, tau=1.0, alpha=0.9,
                  lambda=0.9, Mob=4095) {
  
  ## explain: this function is actually PASL function.

  if (length(beta) != 2) {
    stop("There must be two parameters in the beta vector.")
  }
  f <- beta[1]
  deltaT <- beta[2]

  epoch1 <- TI <= deltaT
  epoch2 <- TI > deltaT & TI < tau + deltaT
  epoch3 <- TI >= tau + deltaT
  
  k <- 1/T1b - 1/T1 - f/lambda
  
  deltaM <- rep(0, length(TI))
  ## fun[epoch1] <- 0
  deltaM[epoch2] <- 2 * Mob * f * (TI[epoch2] - deltaT) * alpha *
    exp(-TI[epoch2] / T1b) * exp(k * TI[epoch2]) *
      (exp(-k * deltaT) - exp(-k * TI[epoch2])) / (k * (TI[epoch2] - deltaT))
  deltaM[epoch3] <- 2 * Mob * f * tau * alpha * exp(-TI[epoch3] / T1b) *
    exp(k * TI[epoch3]) * (exp(-k * deltaT) - exp(-k * (tau + deltaT))) /
      (k * tau) 

  return(deltaM)
  
}

#' @rdname asl2p
#' @export
asl3p <- function(beta, TI, T1b=1.3, T1=1.0, alpha=0.9, lambda=0.9,
                  Mob=100) {

  ## explain: this function is actually PASL function.

  if (length(beta) != 3) {
    stop("There must be three parameters in the beta vector.")
  }
  f <- beta[1]
  deltaT <- beta[2]
  tau <- beta[3]

  ## TI=linspace(0,3.5,701)
  epoch1 <- TI <= deltaT
  epoch2 <- TI > deltaT & TI < (tau + deltaT)
  epoch3 <- TI >= (tau + deltaT)

  k <- 1/T1b - 1/T1 - f/lambda

  fun <- rep(0, length(TI))
  ## fun(epoch1) = 0
  fun[epoch2] <- 2 * Mob * f * (TI[epoch2] - deltaT) * alpha *
    exp(-TI[epoch2] / T1b) * exp(k * TI[epoch2]) *
      (exp(-k * deltaT) - exp(-k * TI[epoch2])) / (k * (TI[epoch2] - deltaT)) 
  fun[epoch3] <- 2 * Mob * f * tau * alpha * exp(-TI[epoch3] / T1b) *
    exp(k * TI[epoch3]) * (exp(-k * deltaT) - exp(-k * (tau + deltaT))) /
      (k * tau)

  return(fun)

}
