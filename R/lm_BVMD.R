
## -------------------------------------------------------------------- ##
#  Function for fitting linear regression model with bivariate response
#  where one of the responses contains missing data.

#  Authors: Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#  Date:    28th of Febryuary 2019
#
#  For details, see "Gangsei, Almøy and Sæbø (2019). Linear Regression 
#  with Bivariate Response Variable Containing Missing Data. 
#  Strategies to Increase Prediction Precision. 

## Arguments ---------------------------------------------------------- ##
#  formula: an object of class "formula" (or one that can be coerced to 
#           that class): a symbolic description of the model to be fitted. 
#  data:    a data frame, list or environment (or object coercible 
#           by as.data.frame to a data frame) containing the variables 
#           in the model. The response variable should be a n1 x 2 
#           matrix, where the last column could/ should contain missing
#           data (NA). Use the 'I()' notation to get a matrix in a data 
#           frame like data <- data.frame(yy = I(),xx1 =,...)
#  method:  The method that should be used for tha analyses 'ML' (maximum
#           likelihood, default) or 'Bayes'
#  hyper:   Only if method is 'Bayes'. A list containing hyperparameters:
#           NULL (Default): Empirical bayes is used to set hyperparameters.
#           If not the following elements are to be implemented:
#           Psi1: a p x p symetric and positive definite matrix
#           Psi2: a p x p symetric and positive definite matrix
#           Phi : a 2 x 2 symetric and positive defenite matrix
#           phi3: a positive scalar (Optional, equal to Phi[1,1] as 
#                  default)
#           zeta1: a scalar > 1.
#           zeta2: a scalar > 2.(Optional, equal to zeta_1 +1 as default)
#           Bet0: a p x 2 matrix of predictor variables


## Value -------------------------------------------------------------- ##
# lm_BVMD returns an object of class "lm_BVMD".
#
# ???? The functions summary and anova are used to obtain and print a 
# summary and analysis of variance table of the results. The generic 
# accessor functions coefficients, effects, fitted.values and residuals 
# extract various useful features of the value returned by lm_BVMD.

# An object of class "lm_BVMD" is a list containing at least the following 
# components:
# 
# formula:        The formula from input
# method:         'ML', 'Bayes', the method used for analysis
# Gram_quad:      A list containg various gram matrices and quadratic sums.
# coefficients:   A list containing coefficient estimates
# residuals:      the residuals, that is response minus fitted values.
# fitted.values:  the fitted mean values.
# FisherInform:   the (2p+3) x (2p+3) observed fisherinformation evaluated
#                 at the "coeffeicients".
# Likelihood:     The log_likelihood when evaluatet at parameter estimates
# JJ33_inv:       The lower right element of inverse fisher information
#                 (only if method is 'ML' or 'FIC').
# om_vec:         Vector used for calculating 'FIC' values in new predictions
#                 (only if method is 'ML' or 'FIC').
# dn_vec:         Vector used for calculating 'FIC' values in new predictions
#                 (only if method is 'ML' or 'FIC').


lm_BVMD <- function(formula,data,method = 'ML',hyper = NULL)
{
## Analysis/ calculations that has to be done for all methods --------- ##
  if(class(formula)!='formula'){formula <- as.formula(formula)}
  
  yy_name <- as.character(lm(formula,data)$terms[[2]])
  data$yy1 <- data[[which(names(data)==yy_name)]][,1]
  data$yy2 <- data[[which(names(data)==yy_name)]][,2]
  
  
  nn1 <- sum(!is.na(data$yy1))
  nn2 <- sum(!is.na(data$yy2))
  ind_n2 <- which(!is.na(data$yy2))
  
  Mod1 <- lm(update.formula(formula,yy1~.),data = data)
  Mod22 <- lm(update.formula(formula,yy2~.),data = data[ind_n2,])
  Mod21 <- lm(update.formula(formula,yy1~.),data = data[ind_n2,])
  
  XX1 <- model.matrix(Mod1)
  XX2 <- model.matrix(Mod22)
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
  qq6 <- 2*QQ[1,1]^2/nn2 + qq5*qq3/nn1
  
## Find coeffecient estimates ----------------------------------------- ##
  # Maximum likelihood ('ML')
  if(method == 'ML')
  {
    Beta1 <- coef(Mod1)
    lambda11 <- nn1/qq3
    lambda22 <- nn2*QQ[1,1]/det(QQ)
    lambda12 <- -QQ[1,2]/QQ[1,1]
    Beta2 <- coef(Mod22) + lambda12*(coef(Mod21)-coef(Mod1)) 
    # Bayes
  }else{
    # Hyper parameters in Empirical Bayes
    if(is.null(hyper))
    {aagg1 <- nlm(Mod_Ev_opt,c(1,1),QQ=as.matrix(qq3),Beta=as.matrix(coef(Mod1)),
                 XtX=XtX1,nn=nn1)
    aagg2 <- nlm(Mod_Ev_opt,c(1,1),QQ=QQ,Beta=
                   as.matrix(cbind(coef(Mod21),coef(Mod22))),XtX=XtX2,nn=nn2)
    
    alpha1 <- aagg1$estimate[1]
    alpha2 <- aagg2$estimate[1]
    gamma1 <- aagg1$estimate[2]
    gamma2 <- aagg2$estimate[2]
    
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
  Likelihood <- pbvlm(theta,XtX1,XtX2,coef(Mod1),
                      cbind(coef(Mod21),coef(Mod22)),QQ,qq3,nn1,nn2)
  
## Observed Fisher information matrix --------------------------------- ##
  FisherInform <- cbind(
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
  
## Inverse observed Fisher information matrix ------------------------- ##
# If method maximum likelihood or FIC, get the lower right element 
#  
#  if((method == 'ML')||(method == 'FIC'))
#  {
#  JJ33_inv <- ((det(QQ)^2/(QQ[1,1]^2*nn2))*
#                 (qq6/(qq6*det(QQ)-2*qq5*QQ[1,2]^2*qq3/nn1
#              -nn2*(qq6*lambda12+(QQ[1,2]/QQ[1,1])*(qq6-(qq5*qq3/nn1)))^2)))
#  }else{
#    JJ33_inv <- NULL
#  }
#
# Get the vector in matrix form used for calculating omegas
#  if((method == 'ML')||(method == 'FIC'))
#  {
#    om_vec <- (diag(rep(1,pp))+
#    as.numeric(QQ[1,1]*nn2^2*qq3*lambda12^2/(det(QQ)*nn1^2))*
#      (nn1/nn2)*XtX2%*%XtX1_inv)%*%as.matrix((coef(Mod1)-coef(Mod21)))
#    Beta2Narr <- as.matrix(coef(Mod22))
    #dn_vec <- sqrt(nn2)*as.matrix(Beta2-Beta2Narr)
#    Dn <- sqrt(nn2)*lambda12
#  }else{
#    om_vec <- NULL
#    Beta2Narr <- NULL
    #dn_vec <- NULL
#    Dn <- NULL
#  }
  

## Predict values ----------------------------------------------------- ## 
  res <- list(formula = formula,
              method = method,
              coefficients = list(
                Beta = cbind(Beta1,Beta2),
                Beta2Narr = as.matrix(coef(Mod22)),
                Beta1Narr = as.matrix(coef(Mod21)),
                Beta1hat = as.matrix(coef(Mod1)),
                Lambda = matrix(c(lambda11,lambda12,lambda12,lambda22),2,2),
                Sigma = matrix(c(1/lambda11,-lambda12/lambda11,
                -lambda12/lambda11,1/lambda22+lambda12^2/lambda11),2,2)
              ),
              Gram_quad = list(QQ=QQ,qq3=qq3,qq4=qq4,qq5=qq5,qq6=qq6,
                               nn1=nn1,nn2=nn2,XtX1 = XtX1,XtX2 = XtX2,
                               XtX1_inv=XtX1_inv,XtX2_inv=XtX2_inv),
              FisherInform = FisherInform,
              Likelihood = Likelihood)#,
           #   JJ33_inv=JJ33_inv,
           #   om_vec = om_vec,
           #   Dn = Dn)
  
  res$model <- Mod1$model
  res$model[[1]] <- data[[which(names(data)==yy_name)]]
  names(res$model)[1] <- yy_name
  res$fitted <- predict.lm_BVMD(res,Newdata = NULL)
  

  # formula:        The formula from input
  # method:         'ML', 'FIC', EB' or 'FB', the method used for analysis
  # coefficients:   A list containing coefficient estimates
  # residuals:      the residuals, that is response minus fitted values.
  # fitted.values:  the fitted mean values.
  # fic.wide:       FIC-value (prediction) for wide model (only if method is 
  #                 'ML' or 'FIC')
  # fic.narrow:     FIC-value for the narrow model (only if method is 'ML' or
  #                 'FIC').
  # FisherInform:   the (2p+3) x (2p+3) observed fisherinformation evaluated
  #                 at the "coeffeicients".
  # Likelihood:     The log_likelihood when evaluatet at parameter estimates
  # JJ33_inv:       The lower right element of inverse fisher information
  #                 (only if method is 'ML' or 'FIC').
  # om_vec:         Vector used for calculating 'FIC' values in new predictions
  #                 (only if method is 'ML' or 'FIC').
  # dn_vec:         Vector used for calculating 'FIC' values in new predictions
  #                 (only if method is 'ML' or 'FIC').
  
  
  
  return(res)
  
}



