
## -------------------------------------------------------------------- ##
#  Function for maximizing model evidxence (Empirical Bayes) in linear 
#  regression model with bivariate response where one of the responses 
#  contains missing data. Used in non linear optimizer "nlm"

#  Authors: Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#  Date:    28th of Febryuary 2019
#
#  For details, see "Gangsei, Almøy and Sæbø (2019). Linear Regression 
#  with Bivariate Response Variable Containing Missing Data. 
#  Strategies to Increase Prediction Precision. 

## Arguments ---------------------------------------------------------- ##
#  aagg:    A vector of length 2 containing the "alpha" and "gamma hyperpar
#  QQ  :    A matrix (2 x 2) or scalar
#  Beta:    Betaparameter
#  XtX:     Gram matrix
#  nn:      Number of observations.
#  ag_max:  Prior set max values for alpha and gamma.

## Value -------------------------------------------------------------- ##
# the negative logarithm of the Model evidence
Mod_Ev_opt <- function(aagg,QQ,Beta,XtX,nn,ag_max=c(Inf,Inf))
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
