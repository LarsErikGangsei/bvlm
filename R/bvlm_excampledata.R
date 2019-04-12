#' Function for generating a data excample used in bvlm.

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
#' @export
bvlm_excampledata <- function()
{

  set.seed(5)
  bvlm_excample <- list(data=NULL,data_test=NULL,Parameters=NULL)
  bvlm_excample$data <- data.frame(Y = I(matrix(0,100,2)),
                     x1 = rnorm(100,2,1),
                     x2 = rgamma(100,3,1),
                     f1 = as.factor(sample(1:3,100,replace=TRUE)))

  set.seed(14)
  bvlm_excample$data_test <- data.frame(Y = I(matrix(0,100,2)),
                                   x1 = rnorm(100,2,1),
                                   x2 = rgamma(100,3,1),
                                   f1 = as.factor(sample(1:3,100,replace=TRUE)))

  set.seed(7)
  bvlm_excample$Parameters$Beta <- matrix(rnorm(10),5,2)%*%matrix(c(1,0.5,0.5,1),2,2)
  colnames(bvlm_excample$Parameters$Beta) <- c('Beta1True','Beta2True')
  Sigma_sq <- matrix(c(4*sqrt(2),1/sqrt(2),1/sqrt(2),2*sqrt(2)),2,2)/sqrt(5)
  Sigma <- Sigma_sq%*%Sigma_sq
  bvlm_excample$Parameters$Lambda <- matrix(c(1/Sigma[1,1],
                                              -Sigma[1,2]/Sigma[1,1],
                                              -Sigma[1,2]/Sigma[1,1],
                                              Sigma[1,1]/det(Sigma)),2,2)
  bvlm_excample$Parameters$Sigma <- Sigma

  set.seed(12)
  E <- matrix(rnorm(200),100,2)%*%Sigma_sq
  bvlm_excample$data$Y <- I(model.matrix(lm(Y~x1+x2+f1,data =
                bvlm_excample$data))%*%bvlm_excample$Parameters$Beta + E)

  set.seed(7)
  E <- matrix(rnorm(200),100,2)%*%Sigma_sq
  bvlm_excample$data_test$Y <- I(model.matrix(lm(Y~x1+x2+f1,data =
                  bvlm_excample$data_test))%*%bvlm_excample$Parameters$Beta + E)


  bvlm_excample$data$Y[21:100,2] <- NA

 return(bvlm_excample)
}

