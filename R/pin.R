#' @title Count error for h2 in breedR
#' 
#' @description 
#' \code{pin} This function counts standard error(se) for heritability(h2) and corr value 
#' and also outputs significent level for corr value in breedR package.
#' 
#' @details 
#' counts standard error(se) for heritability(h2)
#' @aliases pin pin.remlf90
#' @param object an object of breedR result
#' @param formula formula for h2 (or corr, not now)
#' @param signif	 Index to output signif levels, F(default) for non-signif.
#' @param digit a number to control decimal point,5(default).
#' @param vres  Index(T) to return results in vectors, F(default) for direct results.
#' @export sig.level
#' @export
pin <- function(object,formula,signif,digit,vres){
  UseMethod("pin",object)
}

#' @return the result is returned directly.
#' @author Yuanzhen Lin <yzhlinscau@@163.com>
#' @references 
#' breedRPlus website:https://github.com/yzhlinscau/
#' @examples 
#' library(breedR)
#' library(breedRPlus)
#' 
#' data(globulus)
#' res.animal <- remlf90(fixed = phe_X ~ 1,
#'                       random = ~ gg,
#'                       genetic = list(model = 'add_animal',
#'                       pedigree = globulus[, 1:3],
#'                       id = 'self'),
#'                       data = globulus)
#'                       
#' pin(res.animal, h2~V2/(V1+V2+V3))
#' 



#######
# pin.remlf90() function details
######

#' @method  pin remlf90
#' @rdname  pin
#' @export
pin.remlf90 <- function(object,formula,signif=FALSE,digit=5,vres=FALSE) {
  
  #require(msm, quietly = TRUE)
  #require(AAfun)
  
  if (!inherits(object, "breedR")) 
    stop("Argument must be a breedR object")
  
  dd<-gsub('V','x',formula)
  formula<-as.formula(paste(dd[2],dd[3],sep=' ~ '))
  
  transform<-formula
  aa<-object$var[,"Estimated variances"]
  pframe <- as.list(aa)
  names(pframe) <- paste("x", seq(1, length(pframe)), sep = "")
  tvalue<-eval(deriv(transform[[length(transform)]], names(pframe)),pframe)
  tname <- if(length(transform)==3){transform[[2]]}else ""
  
  invAI <- object$reml$invAI
  se <- deltamethod(transform,aa,invAI)
  
  tvalue<-round(tvalue,digit)
  se<-round(se,digit)
  result<-data.frame(row.names=tname, Estimate=tvalue, SE=se)
  result1<-result
  result1$sig.level<-sig.level(tvalue,se)
  
  if(vres==FALSE){cat("\n")
    #options(digits=digit)
    if(signif==TRUE){ 
      print(result1)
      cat("---------------")
      cat("\nSig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n")    
    }else print(result)
    cat("\n")
  }else{
    vv<-vector() #<-NULL
    vv[1]<-tvalue;vv[2]<-se
    names(vv)<-c(tname,paste(tname,"se",sep="."))
    return(vv)
    }
  
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# sig.level functions

sig.level<-function(tvalue,se,...){
  n<-length(se)
  siglevel<-1:n
  for(i in 1:n){    
    sig.se<-c(se[i]*1.450,se[i]*2.326,se[i]*3.090)  
    
    if(abs(tvalue[i])>sig.se[3]) {siglevel[i]<-"***"}
    else if(abs(tvalue[i])>sig.se[2]) {siglevel[i]<-"**"}
    else if(abs(tvalue[i])>sig.se[1]) {siglevel[i]<-"*"}
    else {siglevel[i]<-"Not signif"}
  }
  siglevel
}

##------------------------------
