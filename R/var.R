#' @title Output variance component for BreedR object
#' 
#' @description 
#' \code{var} This function Output variance component for BreedR object.
#' 
#' @details 
#' Output variance component for breedR object,similar to asreml.
#' @aliases var var.remlf90
#' @param object an object of breedR result
#' @param mulT multi-trait model(default, FALSE).

#' @export 
var <- function(object,mulT) {
  UseMethod("var")
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
#' var(res.animal)
#' 


#' @method var remlf90
#' @rdname var
#' @export
var.remlf90 <- function (object,mulT=FALSE) {
  

  
  if (!inherits(object, "breedR")) 
    stop("Argument must be a breedR object")
  
  df<-as.data.frame(summary(object)$var)
  
  df$gamma<-df[,1]/df[nrow(df),1]
  if(mulT==TRUE) df$gamma<-df[,1]
  
  df$z.ratio<-df[,1]/df[,2]
  
  const<-function(x){
    cons.v<-1:length(x)
    for(i in 1:length(x)){
      #if(abs(x[i])!=x[length(x)]) cons.v[i]='Positive'
      if(abs(x[i])<=1e-6) cons.v[i]='Boundary'
      else cons.v[i]='Positive'
    }
    #cons.v[length(x)]='Positive'
    cons.v
  }
  
  df$constraint<-const(df$gamma)
  
  df<-df[,c(3,1:2,4:5)]
  colnames(df)[2:3]<-c('component','std.error')
  
  return(df)
}


