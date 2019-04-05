
## -------------------------------------------------------------------- ##
#  Function returning predicted values for bivariate linear regression 
#  at parameter value theta

#  Authors: Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#  Date:    28th of Febryuary 2019
#
#  For details, see "Gangsei, Almøy and Sæbø (2019). Linear Regression 
#  with Bivariate Response Variable Containing Missing Data. 
#  Strategies to Increase Prediction Precision. 

## Arguments ---------------------------------------------------------- ##
#  Obj:     An object of the lm_BVMD class or a list containing necessarry
#           arguments
#  Newdata: A data frame with new predictors. Names corresponding to the 
#           call/ formula of Obj. 



## Value -------------------------------------------------------------- ##
#  A list containing:
#  fitted:  A matrix of size n x 2 containing fitted values
#  ficvals: If Obj$method is 'ML' or 'FIC', a matrix of size n x 2 
#           containing fic values for narrow and wide model respectively
#  fittedNW:If Obj$method is 'ML' or 'FIC', a matrix of size n x 2
#           containing fitted values for second response using narrow and
#           wide model respectively

predict.lm_BVMD <- function(Obj,Newdata = NULL)
{
## Predictor matrix based on input data and formula
  if(is.null(Newdata))
    {
    Newdata <- as.data.frame(Obj$model)
  }
  
  if(length(as.character(Obj$formula))==3)
  {
    pred_form <- as.formula(paste(as.character(Mod_form)[-2],collapse=''))
  }else
  {
    pred_form <- Obj$formula
  } 

    XXp <- model.matrix(pred_form,data = Newdata)

  fitted <- XXp%*%Obj$coefficients$Beta

## FIC values for prediction
if((Obj$method == 'ML')||(Obj$method == 'FIC'))
{
  #omega_sq <- (XXp%*%Obj$om_vec)^2
  #dn_sq <- (XXp%*%Obj$dn_vec)^2
  #ficvals <- cbind(dn_sq*omega_sq,2*as.numeric(Obj$JJ33_inv)*omega_sq)
  #ficvals <- cbind((Obj$Dn^2)*omega_sq,2*as.numeric(Obj$JJ33_inv)*omega_sq)
  if(Obj$method == 'FIC'){lambda12<-NULL}else{
    (lambda12 <- Obj$coefficients$Lambda[1,2])
  }
  
  ficvals <- apply(XXp,1,lambda12_fit,
                    Beta11h= Obj$coefficients$Beta[,1],
                    Beta21h=Obj$coefficients$Beta1Narr,
                    XtX1_inv = Obj$Gram_quad$XtX1_inv,
                    XtX2 = Obj$Gram_quad$XtX2,
                    QQ=Obj$Gram_quad$QQ,
                    qq3=Obj$Gram_quad$qq3,
                    qq5=Obj$Gram_quad$qq5,
                    qq6=Obj$Gram_quad$qq6,
                    nn1=Obj$Gram_quad$nn1,
                    nn2=Obj$Gram_quad$nn2,
                    lambda12FIC=lambda12)
  ficvals <- matrix(unlist(ficvals),length(ficvals),3,byrow=TRUE)
  
  colnames(ficvals) <- c('lambda12FIC','Narrow','Wide') 
  fittedNW <- cbind(XXp%*%Obj$coefficients$Beta2Narr,fitted[,2])
  colnames(fittedNW) <- c('Narrow','Wide') 
  
  if(Obj$method == 'FIC')
    {
    fittedNW[,2] <- (fittedNW[,1] + ficvals[,1]*
                as.numeric(XXp%*%(Obj$coefficients$Beta1Narr-
                                    Obj$coefficients$Beta[,1])))
    fitted[,2] <- fittedNW[,2]
  fitted[ficvals[,2]<ficvals[,3],2] <- fittedNW[ficvals[,2]<ficvals[,3],1]
  }
}else{
  ficvals <- NULL
  fittedNW <- NULL
}
  
  return(list(fitted = fitted,ficvals = ficvals,fittedNW = fittedNW))
}