#' Function for generating a data excample used in bvlm.
#' @param seed: Used to control "set.seed".
#' @param rhoErr: korrelation in error covariance matrix

#' @return A list bvlm_excample with three elements:
#' @field data: A data frame containing the elements Y(100 x 2), with 85 missing
#' x1 and x2 (both numeric vectors of length 100) and f1 a "factor vector" with
#' 3 levels of length 100
#' @field data_test: Like data, but lacking Y.
#' @field Parameters: A list with 3 elements of "true parameters".
#' Beta (5 x 2), Sigma and Lambda  (both(2 x 2))
#'
#'
#' @author  Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#' @references Gangsei, Almøy and Sæbø (2019). Linear Regression
#'  with Bivariate Response Variable Containing Missing Data.
#'  Strategies to Increase Prediction Precision.
#' @import mvtnorm
#' @export
bvlm_excampledata <- function(seed,rhoErr)
{
  set.seed(seed)
  seeds <- sample(1:100,10)

  set.seed(seeds[1])
  bvlm_excample <- list(data=NULL,data_test=NULL,Parameters=NULL)
  bvlm_excample$data <- data.frame(Y = I(matrix(0,100,2)),
                     x1 = rnorm(100,2,1),
                     x2 = rgamma(100,3,1),
                     f1 = as.factor(sample(1:3,100,replace=TRUE)))


  set.seed(seeds[2])
  bvlm_excample$data_test <- data.frame(Y = I(matrix(0,100,2)),
                                   x1 = rnorm(100,2,1),
                                   x2 = rgamma(100,3,1),
                                   f1 = as.factor(sample(1:3,100,replace=TRUE)))

  set.seed(seeds[3])
  bvlm_excample$Parameters$Beta <- matrix(rnorm(10),5,2)%*%matrix(c(1,0.5,0.5,1),2,2)
  colnames(bvlm_excample$Parameters$Beta) <- c('Beta1True','Beta2True')


  a12 <- ifelse(rhoErr>0,1,-1)*sqrt((1-sqrt(1-rhoErr^2))/2)
  a1 <- sqrt(1-a12^2)


  Sigma_sq <- matrix(c(a1,a12,a12,a1),2,2)
  Sigma <- Sigma_sq%*%Sigma_sq
  bvlm_excample$Parameters$Lambda <- matrix(c(1/Sigma[1,1],
                                              -Sigma[1,2]/Sigma[1,1],
                                              -Sigma[1,2]/Sigma[1,1],
                                              Sigma[1,1]/det(Sigma)),2,2)
  bvlm_excample$Parameters$Sigma <- Sigma
  bvlm_excample$Parameters$Sigma_sq <- Sigma_sq

  set.seed(seeds[4])
  #E <- matrix(rnorm(200),100,2)%*%Sigma_sq
  E <- rmvnorm(100,c(0,0),Sigma)

  bvlm_excample$data$Y <- I(model.matrix(lm(Y~x1+x2+f1,data =
                bvlm_excample$data))%*%bvlm_excample$Parameters$Beta + E)

  set.seed(seeds[9])
  #E <- matrix(rnorm(200),100,2)%*%Sigma_sq
  E <- rmvnorm(100,c(0,0),Sigma)
  bvlm_excample$data_test$Y <- I(model.matrix(lm(Y~x1+x2+f1,data =
                  bvlm_excample$data_test))%*%bvlm_excample$Parameters$Beta + E)


  bvlm_excample$data$Y[21:100,2] <- NA

  colnames(bvlm_excample$data$Y) <- c('y1','y2')
  colnames(bvlm_excample$data_test$Y) <- c('y1','y2')

 return(bvlm_excample)
}

