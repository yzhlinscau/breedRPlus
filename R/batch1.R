#' @title Run batch analysis for single trait in breedR
#' 
#' @description 
#' \code{batch1} This function carries out batch analysis for  
#' mult-trait with same model and also output heritability etc. in breedR package.
#' 
#' @details 
#' breedR batch analysis
#' 
#' @param data   aim dataset
#' @param FMod	 Fixed mode,should be 'y~1+fixed.factors'.
#' @param RMod   Randomed mode, should be '~random.factors'.
#' @param geneticM  Randomed terms for tree model, details see example.
#' @param SpM    Spatial error terms. details see example. 
#' @param pformula formula for h2 (or corr).
#'  
#' @export batch1
#' @return the result is returned directly.
#' @author Yuanzhen Lin <yzhlinscau@@163.com>
#' @references 
#' breedRPlus website:https://github.com/yzhlinscau/
#' @examples 
#' library(breedR)
#' library(breedRPlus)
#' library(tidyr)
#' library(plyr)
#' library(dplyr)
#' 
#' data(douglas)
#' S3<-subset(douglas,site=='s3')
#' #summary(S3);str(S3)
#' 
#' S3a<-dplyr::filter(S3,is.na(dad)) # hs
#' S3a<-transform(S3a,Mum=factor(mum))
#' S3a<-droplevels(S3a)
#' names(S3a)[7:8]<-c('x1','y1')
#' 
#' df<-tidyr::gather(S3a,key=Trait,y,c(-1:-8,-12,-14:-16))
#' #str(df)
#' 
#' 
#' # for parent model
#' fixed = y ~ 1+orig
#' random1=~Mum+block
#' pformula1=h2~4*V1/(V1+V3)
#' 
#' mm<-plyr::ddply(df,'Trait',
#'     function(dat) batch1(dat,FMod=fixed,
#'                              RMod=random1,
#'                              pformula=pformula1)
#'                              )
#' #result                             
#' mm                      
#' 
#' # for tree model
#' random2=~block
#' genetic=list(model = 'add_animal',
#'              pedigree = S3a[,1:3],
#'              id = 'self')
#' pformula2=h2~V2/(V2+V3)
#' 
#' mm1<-plyr::ddply(df,'Trait',
#'     function(dat) batch1(dat,FMod=fixed,
#'                              RMod=random2,
#'                              geneticM=genetic,
#'                              pformula=pformula2)
#'                              )
#' 
#' #result                             
#' mm1
#' 
#'   


# UseMethod("batch1")
batch1 <- function(data,FMod,RMod=NULL,
                       geneticM=NULL,SpM=NULL,
                       pformula=NULL) {
  
  options(warn=-1)
  
  if(!is.null(RMod)){
    if(!is.null(geneticM)){
      if(!is.null(SpM)){
        bdR <- breedR::remlf90(fixed = FMod,
                               random = RMod,
                               genetic=geneticM,
                               spatial=SpM,
                               data = data)
      }else{
        bdR <- breedR::remlf90(fixed = FMod,
                               random = RMod,
                               genetic=geneticM,
                               data = data)
      }
    }else{
      bdR <- breedR::remlf90(fixed = FMod,
                             random = RMod,
                             data = data)
    }
  }else{
    if(!is.null(geneticM)){
      if(!is.null(SpM)){
        bdR <- breedR::remlf90(fixed = FMod,
                               genetic=geneticM,
                               spatial=SpM,
                               data = data)
      }else{
        bdR <- breedR::remlf90(fixed = FMod,
                               genetic=geneticM,
                               data = data)
      }
    }else{
      bdR <- breedR::remlf90(fixed = FMod,
                             data = data)
    }
  }
  
  # get AIC and logLik
  # res.ls is a data frame
  res.ls<-summary(bdR)$model.fit[c(1,3)] 
  
  # get varcomponents, here should be careful
  # transpose var's results (N X 2) into 2 X N
  # N is random factor's total number
  Var<-t(summary(bdR)$var) 
  nn<-ncol(Var) 
  # get var for each random factor
  res.ls[,3:(nn+2)]<-round(Var[1,],3) 
  # get var's se for each random factor
  res.ls[,(nn+3):(2*nn+2)]<-round(Var[2,],3)
  
  NVar<-colnames(Var) # Fam, Residual,etc.
  # rename 'Residual' into 'Rvar'
  #NVar[which(NVar=='Residual')]<-'Rvar' 
  Nvar.se<-paste(NVar,'se',sep='.') # Fam.se, Rvar.se
  names(res.ls)[-1:-2]<-c(NVar,Nvar.se)
  
  #names(vv2) #h2, h2.se
  if(!is.null(pformula)){
    nn1<-ncol(res.ls)
    vv2<-breedRPlus::pin(bdR, pformula,vres=T)   # h2 ~ 4 * V1/(V1+V2)
    res.ls[,nn1+1]<-round(vv2[1],3)   # tvalue
    res.ls[,nn1+2]<-round(vv2[2],3)   # tv.se
    names(res.ls)[-1:-2]<-c(NVar,Nvar.se,names(vv2))
  }
  
  res.ls
}

