##
## demo file for bdR.pin et al. 
##
 library(breedRPlus)
 library(breedR)
 
 data(globulus)
 res.animal <- remlf90(fixed = phe_X ~ 1,
                       random = ~ gg,
                       genetic = list(model = 'add_animal',
                       pedigree = globulus[, 1:3],
                       id = 'self'),
                       data = globulus)
 
 #class(res.animal)
 # output variance components
 var(res.animal)

 # test trait's norm 
 plot1(res.animal)

 # calculate heritability 
 pin(res.animal, h2~V2/(V1+V2+V3))
 
