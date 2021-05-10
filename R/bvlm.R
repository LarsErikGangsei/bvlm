#'  Function for fitting linear regression model with bivariate response
#'  where one of the responses contains missing data.

#' @param formula an object of class "formula" (or one that can be coerced to
#' that class) a symbolic description of the model to be fitted.
#' @param data a data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables
#' in the model. The response variable should be a n1 x 2
#' matrix, where the last column could/ should contain missing
#' data (NA). Use the 'I()' notation to get a matrix in a data
#' frame like data <- data.frame(yy = I(),xx1 =,...)
#' @param method The method that should be used for tha analyses 'ML' (maximum
#' likelihood, default) or 'Bayes'
#' @param hyper Only if method is 'Bayes'. A list containing hyperparameters:
#' NULL (Default): Empirical bayes is used to set hyperparameters.
#' If not the following elements are to be implemented:
#' Psi1: a p x p symetric and positive definite matrix
#' Psi2: a p x p symetric and positive definite matrix
#' Phi : a 2 x 2 symetric and positive defenite matrix
#' phi3: a positive scalar (Optional, equal to Phi[1,1] as default)
#' zeta1: a scalar > 1.
#' zeta2: a scalar > 2.(Optional, equal to zeta_1 +1 as default)
#' Bet0: a p x 2 matrix of predictor variables
#' @param contrasts an optional list. Shoud be set as "contr.sum" for factor variabels
#' if Empirical Bayes is used. See the contrasts.arg of model.matrix.default

#' @return an object of class "bvlm". The generic accessor functions
#' coefficients, fitted.values predict and residuals extract various
#' useful features of the value returned by bvlm.
#'
#' An object of class "bvlm" is a list containing at least the following components:
#'
#' @field formula: the formula used in the function
#' @field contrasts: the contrast used in the function
#' @field method: 'ML', 'Bayes', the method used for analysis
#' @field coefficients: a named list of coefficient estimates
#' @field residuals: the residuals, that is response minus fitted values.
#' @field fitted.values: the fitted mean values.
#' @field gram: A named list containg gram matrices and their inverse for the predictor
#' @field quads: A named list containing different quadratic sums.
#' @field const: A named list containing constants used in FIC prediction
#' @field coefols: Predictor estimates from OLS regression with full and reducet set of data.
#' @field variable matrices based on the full and reduset set of data.
#' @field fisherinform: the (2p+3) x (2p+3) observed fisherinformation evaluated
#' at the coeffeicients.
#' @field likelihood: The log likelihood when evaluatet at parameter estimates
#' @field hyper: Named list of prior parameters. Null if method is 'ML'
#' @field posterior: Named list of posterior parameters. Null if method is 'ML'
#' @author  Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#' @references Gangsei, Almøy and Sæbø (2019). Linear Regression
#'  with Bivariate Response Variable Containing Missing Data.
#'  Strategies to Increase Prediction Precision.
#' @examples
#' # Simulate data
#' bvlm_excample <- bvlm_excampledata(6,0.5)
#'
#' # Fit ML model
#' bvlmObj.ML <- bvlm(formula = 'Y ~ x1 + x2 + f1',
#' data=bvlm_excample$data,method='ML')
#' cbind(bvlm_excample$Parameters$Beta,bvlmObj.ML$coefficients$Beta)[,c(1,3,2,4)]# Compare true and estimated regression parameters
#' cbind(bvlm_excample$Parameters$Sigma,bvlmObj.ML$coefficients$Sigma)#Compare true and estimated error covariance matrix
#' cbind(bvlm_excample$Parameters$Lambda,bvlmObj.ML$coefficients$Lambda)#Compare true and estimated Lambda matrix
#'
#' # Prediction
#' bvlmPrediction.TRUE <- model.matrix(lm('rep(0,100)~ x1 + x2 + f1',data = bvlm_excample$data_test))%*%bvlm_excample$Parameters$Beta# The true expected values
#' bvlmPrediction.ML <- predict_bvlm(bvlmObj.ML,newdata = bvlm_excample$data_test,ficpredict=FALSE)# Prediction on new data, using ML estimates
#' bvlmPrediction.FIC <- predict_bvlm(bvlmObj.ML,newdata = bvlm_excample$data_test,ficpredict=TRUE)# Prediction on new data, using FIC-prediction
#' plot(density(bvlmPrediction.TRUE[,2]-bvlmPrediction.FIC$ficfitted2[,1]),col='red',ylim = c(0,2),
#' main = 'Error densyties. Narrow (red), fic(blue) and ML(green)')
#' points(density(bvlmPrediction.TRUE[,2]-bvlmPrediction.FIC$fitted[,2]),col='blue',type = 'l')
#' points(density(bvlmPrediction.TRUE[,2]-bvlmPrediction.FIC$ficfitted2[,3]),col='green',type = 'l')
#'
#' # RMSEP predicted vs. true prdiction
#' rmsep <- c(sqrt(mean((bvlmPrediction.TRUE[,2]-bvlmPrediction.FIC$ficfitted2[,1])^2)),
#' sqrt(mean((bvlmPrediction.TRUE[,2]-bvlmPrediction.FIC$fitted[,2])^2)),
#' sqrt(mean((bvlmPrediction.TRUE[,2]-bvlmPrediction.ML$fitted[,2])^2)))
#' names(rmsep) <- c('Narrow','FIC','ML')
#' print(rmsep)
#'
#' # Fit Empirical Bayes model
#' bvlmObj.EB <- bvlm(formula = 'Y ~ x1 + x2 + f1',data=bvlm_excample$data,
#' method='Bayes')# Warnings uncentered data and not defined contrast
#' bvlm_excample$data$Y <- scale(bvlm_excample$data$Y,center = TRUE,scale = FALSE)
#' bvlm_excample$data$x1 <- scale(bvlm_excample$data$x1,center = TRUE,scale = TRUE)
#' bvlm_excample$data$x2 <- scale(bvlm_excample$data$x2,center = TRUE,scale = TRUE)
#' bvlm_excample$data_test$x1 <- (bvlm_excample$data_test$x1 - attributes(bvlm_excample$data$x1)$`scaled:center`)/attributes(bvlm_excample$data$x1)$`scaled:scale`
#' bvlm_excample$data_test$x2 <- (bvlm_excample$data_test$x2 - attributes(bvlm_excample$data$x2)$`scaled:center`)/attributes(bvlm_excample$data$x2)$`scaled:scale`
#' bvlmObj.EB <- bvlm(formula = 'Y ~ x1 + x2 + f1',data=bvlm_excample$data,
#' method='Bayes',contrasts = list(f1 = 'contr.sum'))# No warnings
#' bvlmObj.EB$coefficients$Beta[1,] <- bvlmObj.EB$coefficients$Beta[1,] + attributes(bvlm_excample$data$Y)$`scaled:center`
#'
#' # RMSEP predicted vs. true observation
#' bvlmPrediction.EB <- predict_bvlm(bvlmObj.EB,newdata = bvlm_excample$data_test,ficpredict=FALSE)# Prediction on new data, using EB estimates
#' rmsep <- c(sqrt(mean((bvlm_excample$data_test$Y[,2]-bvlmPrediction.TRUE[,2])^2)),
#' sqrt(mean((bvlm_excample$data_test$Y[,2]-bvlmPrediction.FIC$ficfitted2[,1])^2)),
#' sqrt(mean((bvlm_excample$data_test$Y[,2]-bvlmPrediction.FIC$fitted[,2])^2)),
#' sqrt(mean((bvlm_excample$data_test$Y[,2]-bvlmPrediction.FIC$ficfitted2[,3])^2)),
#' sqrt(mean((bvlm_excample$data_test$Y[,2]-bvlmPrediction.EB$fitted[,2])^2)))
#' names(rmsep) <- c('True','Narrow','FIC','ML','EB')
#' print(rmsep)

#' @export
bvlm <- function(formula,data,method = 'ML',hyper = NULL,contrasts=NULL)
  {
  ## Analysis/ calculations that has to be done for all methods --------- ##

  if(class(formula)!='formula'){formula <- as.formula(formula)}
  if(class(data)!='data.frame'){data <- as.data.frame(data)}
  model <- data
  yy_name <- as.character(lm(formula,data)$terms[[2]])
  data$yy1 <- data[[which(names(data)==yy_name)]][,1]
  data$yy2 <- data[[which(names(data)==yy_name)]][,2]

  nn1 <- sum(!is.na(data$yy1))
  nn2 <- sum(!is.na(data$yy2))
  ind_n2 <- which(!is.na(data$yy2))

   #XX2 <- model.matrix(Mod22)

  # If empirical bayes is to be used responses are centered and
  # continious predictors are centered and scaled
  if ((is.null(hyper))&&(method=='Bayes'))
  {
    for(ii in 1:dim(data)[2])
    {
      if (class(data[,ii])=='numeric')
      {
        if (mean(data[,ii])>(10^(-10)))
        {
          warning('The Empirical Bayes strategy is based on centered and scaled (predictors) data')
        }
      }else{
        if (is.null(contrasts))
        {
          warning('The Empirical Bayes strategy is based on centered data and contrasts for factor variables should be set to "contr.sum"')
        }
      }
    }
  }

  Mod1 <- lm(update.formula(formula,yy1~.),data = data,contrasts=contrasts)
  Mod22 <- lm(update.formula(formula,yy2~.),data = data[ind_n2,],contrasts=contrasts)
  Mod21 <- lm(update.formula(formula,yy1~.),data = data[ind_n2,],contrasts=contrasts)

  XX1 <- model.matrix(Mod1)
  XX2 <- XX1[ind_n2,]
  pp <- dim(XX1)[2]

  svdX1 <- svd(XX1)
  svdX2 <- svd(XX2)

  XtX1 <- svdX1$v%*%diag(svdX1$d^2)%*%t(svdX1$v)
  XtX2 <- svdX2$v%*%diag(svdX2$d^2)%*%t(svdX2$v)

  XtX1_inv <- svdX1$v%*%diag(svdX1$d^(-2))%*%t(svdX1$v)
  XtX2_inv <- svdX2$v%*%diag(svdX2$d^(-2))%*%t(svdX2$v)

  QQ <- (rbind(residuals(Mod21),residuals(Mod22))%*%
           cbind(residuals(Mod21),residuals(Mod22)))

  qq3 <- sum(residuals(Mod1)^2)
  qq4 <- (t(as.matrix(coef(Mod1)-coef(Mod21)))%*%XtX2%*%
            as.matrix(coef(Mod1)-coef(Mod21)))
  qq5 <- (t(as.matrix(coef(Mod1)-coef(Mod21)))%*%XtX2%*%XtX1_inv%*%
            XtX2%*%as.matrix(coef(Mod1)-coef(Mod21)))


  ## Find coeffecient estimates ----------------------------------------- ##
  # Maximum likelihood ('ML')
  if(method == 'ML')
  {
    Beta1 <- coef(Mod1)
    lambda11 <- nn1/qq3
    lambda22 <- nn2*QQ[1,1]/det(QQ)
    lambda12 <- -QQ[1,2]/QQ[1,1]
    Beta2 <- coef(Mod22) + lambda12*(coef(Mod21)-coef(Mod1))
    hyper<-NULL
    Upsilon <- NULL
    upsilon3 <- NULL
    # Bayes
  }else{
    # Hyper parameters in Empirical Bayes
    if(is.null(hyper))
    {
      aagg1 <- list(code=5)
      startpos <- c(1,1)
      startiter <- 100
      count <- 0
      options(warn=-1)
      while(aagg1$code>2)
      {
      aagg1 <- nlm(evopt,startpos,QQ=as.matrix(qq3),Beta=as.matrix(coef(Mod1)),
                  XtX=XtX1,nn=nn1,iterlim = startiter)
      if(aagg1$code==5){startpos <- startpos*2}
      if(aagg1$code==3){startpos <- startpos/2}
      if(aagg1$code==4){startiter <- startiter*2}
      count <- count + 1
      if(count>20){break}
      }

      aagg2 <- list(code=5)
      startpos <- c(1,1)
      startiter <- 100
      count <- 0
      while(aagg2$code>2)
      {
        aagg2 <- nlm(evopt,startpos,QQ=QQ,Beta=
                       as.matrix(cbind(coef(Mod21),coef(Mod22))),XtX=XtX2,nn=nn2,
                     iterlim = startiter)
        if(aagg2$code==5){startpos <- startpos*2}
        if(aagg2$code==3){startpos <- startpos/2}
        if(aagg2$code==4){startiter <- startiter*2}
        count <- count + 1
        if(count>20){break()}
      }
      options(warn=0)

    if(aagg1$code<3)
    {
    alpha1 <- aagg1$estimate[1]
    gamma1 <- aagg1$estimate[2]
    }else{
      alpha1 <- 1
      gamma1 <- 1
      warning('Optimalization of model evidence failed, alpa1 and gamma1 set to 1')
    }

    if(aagg2$code<3)
    {
    alpha2 <- aagg2$estimate[1]
    gamma2 <- aagg2$estimate[2]
    }else{
      alpha2 <- 1
      gamma2 <- 1
      warning('Optimalization of model evidence failed, alpa2 and gamma2 set to 1')
    }

  #  print(aagg1$code)
#  print(c(alpha1,gamma1,alpha2,gamma2))
    hyper <- list(Psi1 = (alpha1/nn1)*XtX1,
                  Psi2 = (alpha2/nn2)*XtX2,
                  Phi = diag(rep(gamma2,2)),
                  phi3 = gamma1,
                  zeta1 = gamma1,
                  zeta2 = gamma2,
                  Beta0 = matrix(0,pp,2))
    # print(hyper)
    }
    # Posterior parameters
    Upsilon <- (QQ + hyper$Phi +
                  (rbind(coef(Mod21),coef(Mod22))-t(hyper$Beta0))%*%
                  solve(XtX2_inv+solve(hyper$Psi2))%*%
                  (cbind(coef(Mod21),coef(Mod22))-hyper$Beta0))
    upsilon3 <- (qq3 + hyper$phi3 +
                   t(coef(Mod1)-hyper$Beta0[,1])%*%
                   solve(XtX1_inv+solve(hyper$Psi1))%*%
                   (coef(Mod1)-hyper$Beta0[,1]))

    # Posterior pointestimates (posterior means) Bayesian analysis
    lambda11 <- as.numeric((nn1 +hyper$zeta1)/upsilon3)
    lambda22 <- as.numeric((nn2 + hyper$zeta2)*Upsilon[1,1]/det(Upsilon))
    lambda12 <- -Upsilon[1,2]/Upsilon[1,1]
    Beta1 <- as.vector(solve(hyper$Psi1 + XtX1)%*%
                         (t(XX1)%*%data$yy1+hyper$Psi1%*%hyper$Beta0[,1]))


    names(Beta1) <-names(coef(Mod1))

    Beta2 <- as.vector(solve(hyper$Psi2 + XtX2)%*%
                         (t(XX2)%*%data$yy2[ind_n2]+hyper$Psi2%*%hyper$Beta0[,2])+
                         lambda12*(solve(hyper$Psi2 + XtX2)%*%
                                     (t(XX2)%*%data$yy1[ind_n2]+hyper$Psi2%*%hyper$Beta0[,1])-
                                     Beta1))
    names(Beta2) <-names(coef(Mod22))

  }
  ## Get the likelihood ------------------------------------------------- ##
  theta <- c(lambda11,lambda22,Beta1,Beta2,lambda12)
  likelihood <- pbvlm(theta,XtX1,XtX2,coef(Mod1),
                      cbind(coef(Mod21),coef(Mod22)),QQ,qq3,nn1,nn2)

  ## Observed Fisher information matrix --------------------------------- ##
  fisherinform <- cbind(
    c(nn1/(2*lambda11^2),
      0,
      t(Beta1-coef(Mod1))%*%XtX1,
      rep(0,pp),
      0),
    c(0,
      nn2/(2*lambda22^2),
      rep(0,pp),
      t(Beta2-(coef(Mod22)+lambda12*(coef(Mod21)-coef(Mod1))))%*%XtX2,
      QQ[1,1]*lambda12+QQ[1,2]+
        lambda12*t(Beta1-coef(Mod21))%*%XtX2%*%(Beta1-coef(Mod21))+
        t(Beta2-coef(Mod22))%*%XtX2%*%(Beta1-coef(Mod21))),
    rbind(t(Beta1-coef(Mod1))%*%XtX1,
          rep(0,pp),
          lambda11*XtX1 + lambda12^2*lambda22*XtX2,
          lambda12*lambda22*XtX2,
          lambda22*t(lambda12*(Beta1-coef(Mod21))+Beta2-coef(Mod22))%*%XtX2),
    rbind(rep(0,pp),
          t(Beta2-(coef(Mod22)+lambda12*(coef(Mod21)-coef(Mod1))))%*%XtX2,
          lambda12*lambda22*XtX2,
          lambda22*XtX2,
          lambda22*t(Beta1-coef(Mod21))%*%XtX2),
    c(0,
      QQ[1,1]*lambda12+QQ[1,2]+
        lambda12*t(Beta1-coef(Mod21))%*%XtX2%*%(Beta1-coef(Mod21))+
        t(Beta2-coef(Mod22))%*%XtX2%*%(Beta1-coef(Mod21)),
      lambda22*t(lambda12*(Beta1-coef(Mod21))+Beta2-coef(Mod22))%*%XtX2,
      lambda22*t(Beta1-coef(Mod21))%*%XtX2,
      QQ[1,1]*lambda22+lambda22*t(Beta1-coef(Mod21))%*%XtX2%*%(Beta1-coef(Mod21))
    ))

  ## Return result  ----------------------------------------------------- ##
  res <- list(formula = formula,
              contrasts=contrasts,
              method = method,
              model = model,
              coefficients = list(
                Beta = cbind(Beta1,Beta2),
                Lambda = matrix(c(lambda11,lambda12,lambda12,lambda22),2,2),
                Sigma = matrix(c(1/lambda11,-lambda12/lambda11,
                                 -lambda12/lambda11,1/lambda22+lambda12^2/lambda11),2,2)),
              residuals = cbind(data$yy1,data$yy2)-XX1%*%cbind(Beta1,Beta2),
              fitted.values = XX1%*%cbind(Beta1,Beta2),
              gram = list(XtX1 = XtX1,XtX2 = XtX2,
                          XtX1_inv=XtX1_inv,XtX2_inv=XtX2_inv),
              quads = list(QQ=QQ,qq3=qq3,qq4=qq4,qq5=qq5,nn1=nn1,nn2=nn2),
              coefols  = list(Beta2Narr = as.matrix(coef(Mod22)),
                              Beta1Narr = as.matrix(coef(Mod21)),
                              Beta1Wide = as.matrix(coef(Mod1))),
              fisherinform = fisherinform,
              likelihood = likelihood,
              hyper = hyper,
              posterior = list(Upsilon,upsilon3))#,

  class(res) <- 'bvlm'

  return(res)

}



