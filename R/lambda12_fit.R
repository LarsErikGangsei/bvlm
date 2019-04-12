
#'  Function fitting the value for lambda12 so that FIC value for the
#'  input focus parameter xxN^t x Beta2 is maximized
#'
#' @param xxN: a vector (p x 1) for which the forcus parameter is defined as
#' xxN^t x Beta2. Note that Beta2 is not input in the function.
#' @param Beta11h: OLS estimates for Beta1 based on n1 (all) observations
#' @param Beta21h: OLS estimates for Beta1 based on n2 (full/ not missing data)
#'           observations.
#' @param QQ,qq3,qq5: Quadratic terms, see manuscript for details
#' @param XtX2: Gram matrix for predictors based on full/ not missing data
#' @param XtX1_inv: Inverse gram matrix for predictors based on all data
#' @param nn1,nn2: Total number of observations and number of full observations


#' @return two named vectors:
#' @field lambda12: the lambda12 estimate for the Narrow, FIC and ML situation respectively
#' @field fic: The fic values assosiated with the Narrow, FIC and ML situation.
#' FIC value for narrow is calculated besed on comparisoin with lambda12FIC

#' @author  Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#' @references Gangsei, Almøy and Sæbø (2019). Linear Regression
#'  with Bivariate Response Variable Containing Missing Data.
#'  Strategies to Increase Prediction Precision.

#' @export
lambda12_fit <- function(xxN,Beta11h,Beta21h,XtX1_inv,
                         XtX2,QQ,qq3,qq5,nn1,nn2)
{
  qq6 <- 2*QQ[1,1]^2/nn2 + qq5*qq3/nn1
  DeltaB <- Beta11h-Beta21h
  if(class(xxN)=='matrix'){xxN<-as.numeric(xxN)}

  cc3 <- ((nn2*qq3*QQ[1,1])/(nn1*det(QQ)))*as.numeric(t(as.matrix(xxN))
            %*%XtX2%*%XtX1_inv%*%as.matrix(DeltaB))

  cc4 <- as.numeric(t(as.matrix(xxN))%*%as.matrix(DeltaB))

  cc1 <- det(QQ)-2*QQ[1,2]^2*qq5*qq3/(qq6*nn1)

  cc2 <- -(QQ[1,2]/QQ[1,1])*2*QQ[1,1]^2/(qq6*nn2)

  roots <- polyroot(c(-cc2*cc4,2*cc3*(cc1/(nn2*qq6)-cc2^2)+cc4,
                      3*cc2*cc3,-cc3))

  # polyroot seems to return a small imaginary part even if true
  # real root. Consequently roots <- Re(roots[Im(roots)==0])
  # does not seem to work
  roots <- Re(roots[abs(Im(roots))<(abs(Re(roots))/10^5)])

  if(length(roots)>0)
  {roots <- roots[((-(QQ[1,2]/QQ[1,1]))*roots)>0]}
  if(length(roots)>0)
  {  roots <- roots[abs(roots)<abs(QQ[1,2]/QQ[1,1])]}
  if(length(roots)>1)
  {roots <- roots[abs(roots)==max(abs(roots))]}
  if(length(roots)>0)
  {lambda12FIC <- roots}else
  {lambda12FIC <- -QQ[1,2]/QQ[1,1]}

  FICw <- ((2*det(QQ)^2*(cc3*lambda12FIC^2+cc4)^2)/
        ((cc1-nn2*qq6*(lambda12FIC-cc2)^2)*nn2*QQ[1,1]^2))

  if(FICw<0){FICw<-Inf}

  lambda12ML <- -QQ[1,2]/QQ[1,1]
  FICml <- ((2*det(QQ)^2*(cc3*lambda12ML^2+cc4)^2)/
             ((cc1-nn2*qq6*(lambda12ML-cc2)^2)*nn2*QQ[1,1]^2))

  FICn <- lambda12FIC^2*nn2*(cc3*lambda12FIC^2+cc4)^2

  lambda12_vec <- c(0,lambda12FIC,lambda12ML)
  names(lambda12_vec) <- c('Narrow','FIC','ML')

  fic_vec <- c(FICn,FICw,FICml)
  names(fic_vec) <- c('Narrow','FIC','ML')


  return(list(lambda12=lambda12_vec,
              fic = fic_vec))
}
