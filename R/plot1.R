#' @title Test trait's norm for breedR object
#' 
#' @description 
#' \code{plot1} This function Test trait's norm for breedR object.
#' 
#' @details 
#' Test trait's norm for BreedR object,similar to asreml.
#' @aliases plot1 plot1.remlf90
#' @param object an object of BreedR result
#' @param mulT multi-trait model(default, FALSE).
#' 
#' @export 
plot1 <- function(object,mulT) {
  UseMethod("plot1",object)
}
#' @return the result is returned directly.
#' @author Yuanzhen Lin <yzhlinscau@@163.com>
#' @references 
#' breedRPlus website:https://github.com/yzhlinscau/
#' @examples 
#' library(breedR)
#' library(breedRPlus)
#' 
#' res.animal <- remlf90(fixed = phe_X ~ 1,
#'                       random = ~ gg,
#'                       genetic = list(model = 'add_animal',
#'                       pedigree = globulus[, 1:3],
#'                       id = 'self'),
#'                       data = globulus)
#' plot1(res.animal)
#'



#' @method plot1 remlf90
#' @rdname plot1
#' @export
plot1.remlf90 <- function (object,mulT=FALSE) {
  

  
  if (!inherits(object, "breedR")) 
    stop("Argument must be a breedR object.")
  if(mulT==TRUE)
    stop("Test trait's norm does not works for multi-trait models." )
  
  par(mfrow=c(2,2))
  hist(residuals(object),main='',xlab='Residuals',col='blue')
  qqnorm(residuals(object),main='',col='blue',ylab='Residuals')
  plot(fitted(object),residuals(object),xlab='Fitted',ylab='Residuals',col='blue')
  abline(h=0)
  plot(1:length(fitted(object)),residuals(object),xlab='Unit Number',ylab='Residuals',col='blue')
  abline(h=0)
  par(mfrow=c(1,1))
}

