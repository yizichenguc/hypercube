#####################DATA PREP#######################################
rm(list=ls(all=TRUE))

#install libraries
library(pROC)
library(lmenssp)
library(future.apply)
library(PGEE)

#newdata
setwd("C:\\Users\\CHE4FD\\Documents\\3_BigDataPaper\\code\\yeastG1")
source("lmenssp2.R")
new_solve<-function(A){
  return(chol2inv(chol(A)))}
data(yeastG1)
data <- yeastG1

lmensspRoutine=function(a,significance=NULL){
  #get the correspond data for idx.aux variables, make it into data frame "d1"
  k=length(a)
  colnumber=a+3
  SubsetX=data[,c(1:3,colnumber)]
  
  d1 <- data.frame(SubsetX,SubsetX$time*SubsetX[,-(1:3)])
  nmSubsetX=colnames(SubsetX)
  nmcross<-paste("time_",nmSubsetX[-(1:3)], sep = "")
  colnames(d1)<- c(nmSubsetX,nmcross)
  
  #choose for training data: 80% for training
  idname<-unique(d1$id)
  subject.choose2<-idname[1:226]
  train=d1[d1$id%in%subject.choose2,]
  #test=ctpdata12[ctpdata12$ID%in%subject.choose2,]
  
  ####   second part of the code #####################################
  ###### Model specification ####
  
  #create model formula "temp1"
  fm1<-c("time")
  fm2<-colnames(train[,-(1:3)])
  #fm3<-paste("time_",fm2, sep = "")
  fm<-c(fm1,fm2)
  temp1 <- as.formula(paste("y ~ ", paste(fm, collapse= "+")))
  initial.var<- c(8.273,4.798, 78.837)
  
  #fit the model "lmemssp2", get p-values for all variables
  model1 <- lmenssp2(formula = temp1, data = train,
                     id = train$id, process = "ibm",init=initial.var,
                     timeVar = train$time, silent = TRUE)
  pVals=model1$estimates[3:(3+k-1),4]
  
  #decide which procedure to use for variable selection
  if(significance==2){
    idxSelected=which(pVals %in% sort(pVals)[1:2])
  } else if(significance==1){
    idxSelected=which(pVals %in% sort(pVals)[1])
  }  else{
    idxSelected=which(pVals<significance)
  }
  
  return(list("idxSelected"=idxSelected))
}

#####################################################################################

n = 283
d = 96
seed.HC=NULL

dmHC= 3

#vector.signif=c(2,0.0025,0.01)
#vector.signif=NULL
vector.signif=c(2,0.05)
if(is.null(vector.signif)){
  signif.Default =TRUE
} else{
  signif.Default =FALSE
  #rev:reverse the order of the elements
  vector.signif = rev(vector.signif)
  vector.signif = c(NA,vector.signif)
}

highest.dmHC = dmHC

### Outputs
Matrix.Selection = list()
List.Selection = list()

aux.dmHC4 = aux.dmHC3 = aux.dmHC2 = 'N'

########## case in which dmHC=4 ##########
if(dmHC==4 ){
  
  if(signif.Default==TRUE & highest.dmHC==4){
    aux.signif = 2
  } else if(signif.Default==TRUE & highest.dmHC>4){
    aux.signif = 0.01
  } else{
    aux.signif = vector.signif[4]
  }
  
  if(dmHC==4){
    dimHC = ceiling(d^(1/4))
    nearestHypercube=dimHC^4
    remainderHC=nearestHypercube-d
    if(!is.null(seed.HC)){
      set.seed(seed.HC)
      hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC,dimHC))
    } else{
      hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC,dimHC))
    }}
  
  ### intermediate step to handle 0`s here
  hypercubeSelect=array(0,dim=c(dimHC,dimHC,dimHC,dimHC))
  
  if(all(dim(hypercube)==0)){
    stop(paste('No variables selected at stage 4, increase p-value.'))
  }
  
  for(ind4 in 1:dim(hypercube)[4]){
    for(ind3 in 1:dim(hypercube)[3]){
      for(indR in 1:dim(hypercube)[1]){
        if(length(which(hypercube[indR,,ind3,ind4]!=0))>0){
          index = which(hypercube[indR,,ind3,ind4]!=0)
          idx.aux = hypercube[indR,index,ind3,ind4]
          tryCatch(                       
            expr = {                     
              idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
            },
            
            error = function(e){
              message(paste("There was an error for: (variable id) "),paste(hypercube[indR,,ind3,ind4],collapse = ","))
              idxSelected = NULL
            }
          )
          
          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[indR,idxSelected$idxSelected,ind3,ind4]=hypercubeSelect[indR,idxSelected$idxSelected,ind3,ind4]+1
          }
        }
      } # indR
      for(indC in 1:dim(hypercube)[2]){
        if(length(which(hypercube[,indC,ind3,ind4]!=0))>0){
          index = which(hypercube[,indC,ind3,ind4]!=0)
          idx.aux = hypercube[index,indC,ind3,ind4]
          tryCatch(                       
            expr = {                     
              idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
            },
            
            error = function(e){
              message(paste("There was an error for: (variable id) "),paste(hypercube[,indC,ind3,ind4],collapse = ","))
              idxSelected = NULL
            }
          )
          
          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[idxSelected$idxSelected,indC,ind3,ind4]=hypercubeSelect[idxSelected$idxSelected,indC,ind3,ind4]+1
          }
        }
      } # indC
    } # ind3
    for(indR in 1:dim(hypercube)[1]){
      for(indC in 1:dim(hypercube)[2]){
        if(length(which(hypercube[indR,indC,,ind4]!=0))>0){
          index = which(hypercube[indR,indC,,ind4]!=0)
          idx.aux = hypercube[indR,indC,index,ind4]
          tryCatch(                       
            expr = {                     
              idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
            },
            
            error = function(e){
              message(paste("There was an error for: (variable id) "),paste(hypercube[indR,indC,,ind4],collapse = ","))
              idxSelected = NULL
            }
          )
          
          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[indR,indC,idxSelected$idxSelected,ind4]=hypercubeSelect[indR,indC,idxSelected$idxSelected,ind4]+1
          }
        }
      } # indC
    } # indR
  } # ind4
  # now traverse in the 4th dimension
  for(ind3 in 1:dim(hypercube)[3]){
    for(indR in 1:dim(hypercube)[1]){
      for(indC in 1:dim(hypercube)[2]){
        if(length(which(hypercube[indR,indC,ind3,]!=0))>0){
          index = which(hypercube[indR,indC,ind3,]!=0)
          idx.aux = hypercube[indR,indC,ind3,index]
          tryCatch(                       
            expr = {                     
              idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
            },
            
            error = function(e){
              message(paste("There was an error for: (variable id) "),paste(hypercube[indR,indC,ind3,],collapse = ","))
              idxSelected = NULL
            }
          )
          
          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[indR,indC,ind3,idxSelected$idxSelected]=hypercubeSelect[indR,indC,ind3,idxSelected$idxSelected]+1
          }
        }
      } # indC
    } # indR
  } # ind3
  
  setSelected4Times = hypercube[which(hypercubeSelect>3, arr.ind = TRUE)]
  setSelected3Times = hypercube[which(hypercubeSelect>2, arr.ind = TRUE)]
  setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
  setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]
  
  numSelected4times5=length(setSelected4Times)
  numSelected3times5=length(setSelected3Times)
  numSelected2times5=length(setSelected2Times)
  numSelected1times5=length(setSelected1Times)
  
  Matrix.Selection[[paste('Hypercube with dim',4)]] = c(numSelected1times5,numSelected2times5,numSelected3times5,numSelected4times5)
  names(Matrix.Selection[[paste('Hypercube with dim',4)]]) = c('numSelected1','numSelected2','numSelected3','numSelected4')
  
  List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected1')]] = setSelected1Times
  List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected2')]] = setSelected2Times
  List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected3')]] = setSelected3Times
  List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected4')]] = setSelected4Times
  
  aux.dmHC4 <- 'Y' #readline(cat("Reduction of dimension 4 done!", "\n", length(setSelected4Times5),
  #        "Variables selected at least 4 times","\n",length(setSelected3Times5),
  #        "Variables selected at least 3 times","\n",length(setSelected2Times5),
  #       "Variables selected at least 2 times","\n", "Wanna proceed with reduction?[Y/N]"))
}

########## case in which dmHC=3 ##########
if(dmHC==3 | aux.dmHC4=='Y'){
  
  if(signif.Default==TRUE & highest.dmHC==3){
    aux.signif = 2
  } else if(signif.Default==TRUE & highest.dmHC>3){
    aux.signif = 0.01
  } else{
    aux.signif = vector.signif[3]
  }
  
  if(dmHC==3){
    dimHC = ceiling(d^(1/3))
    nearestHypercube=dimHC^3
    remainderHC=nearestHypercube-d
    if(!is.null(seed.HC)){
      set.seed(seed.HC)
      #arrange variable id in hypercube in a random arrange manner
      hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC))
    } else{
      hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC))
    }} else if(dmHC>3){
      dimHC = ceiling(length(setSelected2Times)^(1/3))
      nearestHypercube=dimHC^3
      remainderHC=nearestHypercube-length(setSelected2Times)
      if(!is.null(seed.HC)){
        set.seed(seed.HC)
        hypercube=array(sample(c(setSelected2Times,rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC))
      } else{
        hypercube=array(c(setSelected2Times,rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC))
      }
    }
  
  ### intermediate step to handle 0`s here
  hypercubeSelect=array(0,dim=c(dimHC,dimHC,dimHC))
  
  if(all(dim(hypercube)==0)){
    stop(paste('No variables selected at stage 3, increase p-value.'))
  }
  
  for(indL in 1:dim(hypercube)[3]){
    for(indR in 1:dim(hypercube)[1]){
      if(length(which(hypercube[indR,,indL]!=0))>0){
        index = which(hypercube[indR,,indL]!=0)
        idx.aux = hypercube[indR,index,indL]
        tryCatch(                       
          expr = {                     
            idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
          },
          
          error = function(e){
            message(paste("There was an error for: (variable id) "),paste(hypercube[indR,,indL],collapse = ","))
            idxSelected = NULL
          }
        )
        
        if(length(idxSelected$idxSelected)>0){
          hypercubeSelect[indR,idxSelected$idxSelected,indL]=hypercubeSelect[indR,idxSelected$idxSelected,indL]+1
        }
      }
    } # indR
    for(indC in 1:dim(hypercube)[2]){
      if(length(which(hypercube[,indC,indL]!=0))>0){
        index = which(hypercube[,indC,indL]!=0)
        idx.aux = hypercube[index,indC,indL]
        tryCatch(                       
          expr = {                     
            idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
          },
          
          error = function(e){
            message(paste("There was an error for: (variable id) "),paste(hypercube[,indC,indL],collapse = ","))
            idxSelected = NULL
          }
        )
        
        if(length(idxSelected$idxSelected)>0){
          hypercubeSelect[idxSelected$idxSelected,indC,indL]=hypercubeSelect[idxSelected$idxSelected,indC,indL]+1
        }
      }
    } # indC
  } # indL
  for(indR in 1:dim(hypercube)[1]){
    for(indC in 1:dim(hypercube)[2]){
      if(length(which(hypercube[indR,indC,]!=0))>0){
        index = which(hypercube[indR,indC,]!=0)
        idx.aux = hypercube[indR,indC,index]
        tryCatch(                       
          expr = {                     
            idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
          },
          
          error = function(e){
            message(paste("There was an error for: (variable id) "),paste(hypercube[indR,indC,],collapse = ","))
            idxSelected = NULL
          }
        )
        
        if(length(idxSelected$idxSelected)>0){
          hypercubeSelect[indR,indC,idxSelected$idxSelected]=hypercubeSelect[indR,indC,idxSelected$idxSelected]+1
        }
      }
    } # indC
  } # indR
  
  setSelected3Times = hypercube[which(hypercubeSelect>2, arr.ind = TRUE)]
  setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
  setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]
  
  numSelected3times=length(setSelected3Times)
  numSelected2times=length(setSelected2Times)
  numSelected1times=length(setSelected1Times)
  
  Matrix.Selection[[paste('Hypercube with dim',3)]] = c(numSelected1times,numSelected2times,numSelected3times)
  names(Matrix.Selection[[paste('Hypercube with dim',3)]]) = c('numSelected1','numSelected2','numSelected3')
  
  if(numSelected1times>1){
    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected1')]] = setSelected1Times
  } else{
    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected1')]] = c(0,ifelse(length(setSelected1Times)==0,0,c(setSelected1Times)))
    List.Selection$`Hypercube with dim 3`$numSelected1 = List.Selection$`Hypercube with dim 3`$numSelected1[-1]
  }
  
  if(numSelected2times>1){
    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected2')]] = setSelected2Times
  } else{
    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected2')]] = c(0,ifelse(length(setSelected2Times)==0,0,c(setSelected2Times)))
    List.Selection$`Hypercube with dim 3`$numSelected2 = List.Selection$`Hypercube with dim 3`$numSelected2[-1]
  }
  
  if(numSelected3times>1){
    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected3')]] = setSelected3Times
  } else{
    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected3')]] = c(0,ifelse(length(setSelected3Times)==0,0,c(setSelected3Times)))
    List.Selection$`Hypercube with dim 3`$numSelected3 = List.Selection$`Hypercube with dim 3`$numSelected3[-1]
  }
  
  aux.dmHC3 <- 'Y'
  
}

########## case in which dmHC=2 ##########
if(dmHC==2 | aux.dmHC3=='Y'){
  
  if(signif.Default==TRUE & highest.dmHC==2){
    aux.signif = 2
  } else if(signif.Default==TRUE & highest.dmHC>2){
    aux.signif = 0.01
  } else{
    aux.signif = vector.signif[2]
  }
  
  if(dmHC==2){
    dimHC = ceiling(d^(1/2))
    nearestHypercube=dimHC^2
    remainderHC=nearestHypercube-d
    if(!is.null(seed.HC)){
      set.seed(seed.HC)
      hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC))
    } else{
      hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC))
    }} else if(dmHC>2){
      dimHC = ceiling(numSelected1times^(1/2))
      nearestHypercube=dimHC^2
      remainderHC=nearestHypercube-length(setSelected2Times)
      if(!is.null(seed.HC)){
        set.seed(seed.HC)
        hypercube=array(sample(c(setSelected2Times,rep(0,remainderHC))),dim=c(dimHC,dimHC))
      } else{
        hypercube=array(c(setSelected2Times,rep(0,remainderHC)),dim=c(dimHC,dimHC))
      }
    }
  
  ### intermediate step to handle 0`s here
  hypercubeSelect=array(0,dim=c(dimHC,dimHC))
  
  if(all(dim(hypercube)==0)){
    stop(paste('No variables selected at stage 2, increase p-value.'))
  }
  
  for(indR in 1:dim(hypercube)[1]){
    if(length(which(hypercube[indR,]!=0))>0){
      index = which(hypercube[indR,]!=0)
      idx.aux = hypercube[indR,index]
      tryCatch(                       
        expr = {                     
          idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
        },
        
        error = function(e){
          message(paste("There was an error for: (variable id) "),paste(hypercube[indR,],collapse = ","))
          idxSelected = NULL
        }
      )
      
      if(length(idxSelected$idxSelected)>0){
        hypercubeSelect[indR,idxSelected$idxSelected]=hypercubeSelect[indR,idxSelected$idxSelected]+1
      }
    }
  } # indR
  for(indC in 1:dim(hypercube)[2]){
    if(length(which(hypercube[,indC]!=0))>0){
      index = which(hypercube[,indC]!=0)
      idx.aux = hypercube[index,indC]
      tryCatch(                       
        expr = {                     
          idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
        },
        
        error = function(e){
          message(paste("There was an error for: (variable id) "),paste(hypercube[,indC],collapse = ","))
          idxSelected = NULL
        }
      )
      
      if(length(idxSelected$idxSelected)>0){
        hypercubeSelect[idxSelected$idxSelected,indC]=hypercubeSelect[idxSelected$idxSelected,indC]+1
      }
    }
  } # indC
  
  setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
  setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]
  
  numSelected2times=length(setSelected2Times)
  numSelected1times=length(setSelected1Times)
  
  Matrix.Selection[[paste('Hypercube with dim',2)]] = c(numSelected1times,numSelected2times)
  names(Matrix.Selection[[paste('Hypercube with dim',2)]]) = c('numSelected1','numSelected2')
  
  if(numSelected1times>1){
    List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected1')]] = setSelected1Times
  } else{
    List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected1')]] = c(0,ifelse(length(setSelected1Times)==0,0,c(setSelected1Times)))
    List.Selection$`Hypercube with dim 2`$numSelected1 = List.Selection$`Hypercube with dim 2`$numSelected1[-1]
  }
  
  if(numSelected2times>1){
    List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected2')]] = setSelected2Times
  } else{
    List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected2')]] = c(0,ifelse(length(setSelected2Times)==0,0,c(setSelected2Times)))
    List.Selection$`Hypercube with dim 2`$numSelected2 = List.Selection$`Hypercube with dim 2`$numSelected2[-1]
  }
  
  aux.dmHC2 <- 'Y' #readline(cat("Reduction of dimension 2 done!", "\n", length(setSelected2Times),
  
}

list("Matrix.Selection" = Matrix.Selection, "List.Selection" = List.Selection)

save(Matrix.Selection, file = 'Matrix.Selection.RData')
save(List.Selection, file = 'List.Selection.RData')