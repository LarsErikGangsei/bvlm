
## -------------------------------------------------------------------- ##
#  Function fitting the value for lambda12 so that FIC value for the 
#  input focus parameter xxN^t x Beta2 is maximized

#  Authors: Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#  Date:    28th of Febryuary 2019
#
#  For details, see "Gangsei, Almøy and Sæbø (2019). Linear Regression 
#  with Bivariate Response Variable Containing Missing Data. 
#  Strategies to Increase Prediction Precision. 

## Arguments ---------------------------------------------------------- ##
#  xxN: a vector (p x 1) for which the forcus parameter is defined as 
#       xxN^t x Beta2. Note that Beta2 is not input in the function.
#  Beta11h: OLS estimates for Beta1 based on n1 (all) observations
#  Beta21h: OLS estimates for Beta1 based on n2 (full/ not missing data) 
#           observations.  
#  QQ, qq3,qq5,qq6: Quadratic terms, see manuscript for details
#  XtX2: Gram matrix for predictors based on full/ not missing data 
#  XtX1_inv: Inverse gram matrix for predictors based on all data 
#  nn1,nn2: Total number of observations and number of full observations


## Value -------------------------------------------------------------- ##
#  A list containing:
#  lambda12FIC:  The fitted lambda12 value
#  FICw: The FIC value assosiated with the wide model using lambda12FIC

lambda12_fit <- function(xxN,Beta11h,Beta21h,XtX1_inv,
                         XtX2,QQ,qq3,qq5,qq6,nn1,nn2,lambda12FIC=NULL)
{
  DeltaB <- Beta11h-Beta21h
  if(class(xxN)=='matrix'){xxN<-as.numeric(xxN)}
  
  cc3 <- ((nn2*qq3*QQ[1,1])/(nn1*det(QQ)))*as.numeric(t(as.matrix(xxN))
            %*%XtX2%*%XtX1_inv%*%as.matrix(DeltaB))
  
  cc4 <- as.numeric(t(as.matrix(xxN))%*%as.matrix(DeltaB))
  
  cc1 <- det(QQ)-2*QQ[1,2]^2*qq5*qq3/(qq6*nn1)
  
  cc2 <- -(QQ[1,2]/QQ[1,1])*2*QQ[1,1]^2/(qq6*nn2)
  
  roots <- polyroot(c(-cc2*cc4,2*cc3*(cc1/(nn2*qq6)-cc2^2)+cc4,
                      3*cc2*cc3,-cc3))
  
  #print(abs(Im(roots)))
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
  {lambda12FIC <- 0}
  
  
  FICw <- ((2*det(QQ)^2*(cc3*lambda12FIC^2+cc4)^2)/
        ((cc1-nn2*qq6*(lambda12FIC-cc2)^2)*nn2*QQ[1,1]^2))
  
  if(FICw<0){FICw<-Inf}
  
  FICn <- lambda12FIC^2*nn2*(cc3*lambda12FIC^2+cc4)^2
  
  return(list(lambda12FIC=as.numeric(lambda12FIC),
              FICn=as.numeric(FICn),
              FICw=as.numeric(FICw)))
}