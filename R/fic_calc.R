
#' Function returning FIC values for wide and narrow model using ML
#' estimates in the bivariate regression model exept for lambda12, which
#' is an argument to the function. Used for making figures.

#' @param lambda12: The value for lambda12 for which to evaluate the FIC values.
#' @param xxN: a vector (p x 1) for which the forcus parameter is defined as
#'       xxN^t x Beta2. Note that Beta2 is not input in the function.
#' @param Beta11h: OLS estimates for Beta1 based on n1 (all) observations
#' @param Beta21h: OLS estimates for Beta1 based on n2 (full/ not missing data)
#'           observations.
#' @param QQ, qq3,qq5,qq6: Quadratic terms, see manuscript for details
#' @param XtX2: Gram matrix for predictors based on full/ not missing data
#' @param XtX1_inv: Inverse gram matrix for predictors based on all data
#' @param nn1,nn2: Total number of observations and number of full observations


#' @return A vector of length two containing two elements, the FIC vales for the
#  wide and narrow model respectively.

#' @author  Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#' @references Gangsei, Almøy and Sæbø (2019). Linear Regression
#'  with Bivariate Response Variable Containing Missing Data.
#'  Strategies to Increase Prediction Precision.

#' @export

fic_calc <- function(lambda12,xxN,Beta11h,Beta21h,XtX1_inv,
                         XtX2,QQ,qq3,qq5,qq6,nn1,nn2)
{
  pp <- length(xxN)

  omega <- (-t((diag(rep(1,pp))+
               as.numeric(QQ[1,1]*nn2^2*qq3*lambda12^2/(det(QQ)*nn1^2))*
               (nn1/nn2)*XtX2%*%XtX1_inv)%*%as.matrix(Beta21h-Beta11h))
                %*%as.matrix(xxN))

  cc1 <- det(QQ)-2*QQ[1,2]^2*qq5*qq3/(qq6*nn1)

  cc2 <- -(QQ[1,2]/QQ[1,1])*2*QQ[1,1]^2*nn1/(2*QQ[1,1]^2*nn1 + nn2*qq5*qq3)

  jj1_inv <- (det(QQ)^2/((cc1-nn2*qq6*(lambda12-cc2)^2)*nn2*QQ[1,1]^2))


  FICw <- 2*omega^2*jj1_inv
  FICn <- lambda12^2*nn2*omega^2

  result <- c(FICw,FICn)
  names(result) <- c('FICw','FICn')

  return(result)
}
