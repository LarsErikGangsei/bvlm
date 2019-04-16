
#'  Function returning predicted values for bivariate linear regression

#' @param Obj: An object of class 'bvlm'
#' @param newdata: A data frame containing observations that shall be predicted.
#'  Names corresponding to the formula of Obj.
#' @param ficpredict: Logical value indicating if fic prediction shall be applied.
#' Default is FALSE. Only possible if Obj$method = 'ML'

#' @return A named list
#' @field fitted: a n x 2 matrix with fitted values for both responses
#' @field lambda12:  Only if ficpredict is TRUE. A n x 3 matrix fic
#' values for each new prediction assosiated with  Narrow, FIC and ML situation.
#' @field fic: Only if ficpredict is TRUE. A n x 3 matrix containing
#' fic values for each new prediction assosiated with  Narrow, FIC
#' and ML situation.
#' @field ficfitted2: Only if ficpredict is TRUE. A n x 3 matrix containing
#' predicted values for the second response assosiated with the  Narrow, FIC
#' and ML situation.
#'
#' @author  Lars Erik Gangsei, Trygve Almøy and Solve Sæbø
#' @references Gangsei, Almøy and Sæbø (2019). Linear Regression
#'  with Bivariate Response Variable Containing Missing Data.
#'  Strategies to Increase Prediction Precision.

#' @export
predict_bvlm <- function(Obj,newdata = NULL,ficpredict=FALSE)
{
## Predictor matrix based on input data and formula
  if(is.null(newdata))
    {
    newdata <- as.data.frame(Obj$model)
  }

  if(length(as.character(Obj$formula))==3)
  {
    pred_form <- as.formula(paste(as.character(Obj$formula)[-2],collapse=''))
  }else
  {
    pred_form <- Obj$formula
  }

    XXp <- model.matrix(pred_form,contrasts = Obj$contrasts,data = newdata)


  fitted <- XXp%*%Obj$coefficients$Beta

  colnames(fitted) <- c('y1','y2')

  res <- list(fitted = fitted)

## FIC values for prediction
if((Obj$method == 'ML')&&(ficpredict == TRUE))
{

  ficvals <- apply(XXp,1,lambda12_fit,
                    Beta11h= Obj$coefols$Beta1Wide,
                    Beta21h=Obj$coefols$Beta1Narr,
                    XtX1_inv = Obj$gram$XtX1_inv,
                    XtX2 = Obj$gram$XtX2,
                    QQ=Obj$quads$QQ,
                    qq3=Obj$quads$qq3,
                    qq5=Obj$quads$qq5,
                    nn1=Obj$quads$nn1,
                    nn2=Obj$quads$nn2)

   ficvals <- matrix(unlist(ficvals),length(ficvals),6,byrow=TRUE)

  lambda12 <- ficvals[,1:3]
  colnames(lambda12) <- c('Narrow','FIC','ML')

  fic <- ficvals[,4:6]
  colnames(fic) <- c('Narrow','FIC','ML')

  fic_extra <- XXp%*%(Obj$coefols$Beta1Narr-Obj$coefols$Beta1Wide)

  ficfitted2 <- cbind(XXp%*%Obj$coefols$Beta2Narr,
                      XXp%*%Obj$coefols$Beta2Narr+lambda12[,2]*fic_extra,
                      XXp%*%Obj$coefols$Beta2Narr+lambda12[,3]*fic_extra)

  i_mat <- as.matrix(apply(fic,1,min))%*%t(as.matrix(rep(1,3)))==fic
  i_mat <- i_mat/(as.matrix(rowSums(i_mat))%*%t(as.matrix(rep(1,3))))

  fitted[,2] <-rowSums(ficfitted2*(i_mat))

  res <- list(fitted = fitted,lambda12 = lambda12,
              fic = fic,ficfitted2 = ficfitted2)
}

  return(res)
}
