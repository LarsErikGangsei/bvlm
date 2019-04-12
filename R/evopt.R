
#'  Function for maximizing model evidxence (Empirical Bayes) in linear
#' regression model with bivariate response where one of the responses
#' contains missing data. Used as input in non linear optimizer "nlm"

#' @param aagg    A vector of length 2 containing the "alpha" and "gamma hyperpar
#' @param QQ      A matrix (2 x 2) or scalar
#' @param Beta    Betaparameter
#' @param XtX     Gram matrix
#' @param nn      Number of observations.
#' @param ag_max  Prior set max values for alpha and gamma.

#' @return the negative logarithm of the Model evidence avaluated at with
#' hyperparameters alpha and gamma given in aagg

#' @author  Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#' @references Gangsei, Almøy and Sæbø (2019). Linear Regression
#'  with Bivariate Response Variable Containing Missing Data.
#'  Strategies to Increase Prediction Precision.

#' @export
evopt <- function(aagg,QQ,Beta,XtX,nn,ag_max=c(Inf,Inf))
{#Start function
  alpha <- aagg[1]
  gamma <- aagg[2]
  pp <- dim(XtX)[1]
  Phi <- diag(gamma,dim(QQ)[1])
  Ups <- QQ + Phi + t(Beta)%*%((alpha/(nn+alpha))*XtX)%*%Beta
  if ((alpha>0)&&(gamma>0)&&
      (alpha<ag_max[1])&&(gamma<ag_max[2]))
  {res <- (lgamma((nn+gamma)/2)
           -lgamma(gamma/2)
           +((gamma+1-dim(QQ)[1])/2)*log(gamma)
           -((nn + gamma -1 +dim(QQ)[1])/2)*log(det(Ups))
           +ifelse(dim(QQ)[1]==2,((nn+gamma-1)/2)*log(Ups[1,1]),0)
           +(pp/2)*(log(alpha)-log(nn+alpha)))
          }else{
             res <- -Inf}
  return(-res)
}#End function
