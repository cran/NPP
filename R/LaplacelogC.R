

LaplacelogC <- function(delta, loglikmle, detHessian, ntheta){
  loglik2 <- delta*loglikmle
  adj <- (ntheta/2)*log(2*pi)-delta*detHessian/2
  return(loglik2+adj)
}
