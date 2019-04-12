#'  Function returning the probability density value (pdf-value) for bivariate linear regression at parameter value theta

#' @param theta:   A vector of length 2p + 3 containing values for lambda11,
#           lambda22,Beta1,Beta2 and lambda12
#' @param XtX1:    Gram matrix (p x p) based on X1.
#' @param XtX2:    Gram matrix (p x p) based on X2.
#' @param Betah1:  (p x 1) ols estimate for Beta1
#' @param Betah2:  (px 2) ols estimate for Beta1 and Beta2 based on n2
#' @param QQ:      2 x 2 squared residual matrix
#' @param qq3:     SSE for linear model based on XX1
#' @param nn1:     Number of observations (X1).
#' @param nn2:     Number of observations (X2).
#' @param log_r:   Logical. When TRUE (default) log likelihood is returned.

#' @return The likelihood (log likelihood is default) for model
#' and data in question avaluated at parameterestimates defined in theta

#' @author  Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#' @references Gangsei, Almøy and Sæbø (2019). Linear Regression
#'  with Bivariate Response Variable Containing Missing Data.
#'  Strategies to Increase Prediction Precision.

#' @export

pbvlm <- function(theta,XtX1,XtX2,Betah1,Betah2,QQ,qq3,nn1,nn2,log_r = TRUE)
{#Start function
  pp <- dim(XtX1)[1]
  lambda11 <- theta[1]
  lambda22 <- theta[2]
  Beta1 <- theta[3:(2+pp)]
  Beta2 <- theta[(3+pp):(2*pp+2)]
  lambda12 <- theta[2*pp+3]

  Omega_inv <- cbind(rbind(lambda11*XtX1+lambda12^2*lambda22*XtX2,
                           lambda12*lambda22*XtX2),
                     rbind(lambda12*lambda22*XtX2,
                           lambda22*XtX2))

  Delta_B <- as.matrix(c(Beta1-Betah1,
                         Beta2-Betah2[,2]-lambda12*(Betah2[,1]-Betah1)))


  res <- (0.5*(-(nn1+nn2)*(log(2)+log(pi)) +
          nn1*log(lambda11) +
            nn2*log(lambda22) -
            qq3*lambda11 -
            QQ[1,1]*lambda22*lambda12^2 -
            2*QQ[1,2]*lambda12*lambda22-
            QQ[2,2]*lambda22 -
            t(Delta_B)%*%Omega_inv%*%Delta_B))

  if(log_r==FALSE)
  {
    res <- exp(res)
  }

  return(res)
}#End function
