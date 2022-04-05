#####################DATA PREP#######################################
rm(list=ls(all=TRUE))

#install libraries
library(pROC)
library(lmenssp)
library(future.apply)
library(ggplot2) 
library(dplyr) 

#set working directory to user defined path
setwd("usr\\bin\\hypercube")

#load simulated data
load("data_sim")
result <- data_sim

#function for Choleski factorization
new_solve<-function(A){
  return(chol2inv(chol(A)))}

#key longitudinal model fitting function
lmenssp2 <-
  function(formula, data = NULL, id, process = "bm", timeVar, init = NULL, tol = 1e-5, maxiter = 100, silent = TRUE){ 
    
    #library(nlme)
    #library(MASS)
    #library(geoR)
    
    mf <- model.frame(formula = formula, data = data)
    y  <- as.matrix(model.extract(mf, "response"))
    x  <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))
    colnames(x)[1] <- gsub("[[:punct:]]", "", colnames(x)[1])
    
    nsubj  <- length(unique(id)) # number of subjects
    ntotal <- nrow(y)            # total number of observations
    #idlist <- unique(id)
    
    Time  <- tapply(timeVar, id, function(x) x)
    Intercept <- rep(1, nrow(data))
    data2 <- data.frame(cbind(id, x))
    DM    <- split(data2[, -1], data2$id)
    YM    <- tapply(y, id, function(x) x)
    nobs  <- tapply(id, id, function(x) length(x)) # number of observations per patient
    
    
    ##########################################################
    ################### INTEGRATED BROWNIAN MOTION ###########
    ##########################################################
    
    if(process == "ibm"){
      
      if (length(init) == 0){
        data.init         <- data
        data.init$timeVar <- timeVar
        data.init$id      <- id
        init <- as.numeric(VarCorr(lme(formula, random = ~ timeVar|id, method = "ML", data = data.init))[,1])
      }
      
      theta.new <- init
      theta.old <- init * 20
      tol       <- tol
      iter      <- 1
      Niter     <- maxiter
      
      
      while (sqrt((theta.old - theta.new) %*% (theta.old - theta.new)) > tol & iter <= Niter){
        
        theta.old <- theta.new
        
        a <- theta.old[1]  ## omegasq
        b <- theta.old[2]  ## sigmasq
        c <- theta.old[3]  ## nu
        
        
        ### betahat
        
        sum.left.beta <- sum.right.beta <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv  <- new_solve(Vi)
          Xi.transp <- t(Xi)
          
          sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
          sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi
          
        }#for (i in 1 : nsubj)
        
        b.hat <- new_solve(sum.left.beta) %*% sum.right.beta
        
        ####### score for theta
        
        ### a = omegasq
        
        sum.a <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          ri       <- Yi - Xi %*% b.hat 
          
          derVa <- Ji
          
          sum.a <- sum.a + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVa))
          
        }#for (i in 1 : nsubj)
        
        
        ### b = sigmasq
        
        sum.b <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          ri       <- Yi - Xi %*% b.hat 
          
          derVb <- Ri
          
          sum.b <- sum.b + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVb))
          
        }#for (i in 1 : nsubj)
        
        
        # c = tausq
        
        sum.c <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          ri       <- Yi - Xi %*% b.hat 
          
          derVc <- Ii
          
          sum.c  <- sum.c + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVc))
          
        }#for (i in 1 : nsubj)
        
        score.theta <- 0.5 * matrix(c(sum.a, sum.b, sum.c))
        
        
        ######### INFORMATION MATRIX
        
        #### 1) a = omegasq, 2) a = omegasq
        
        sum.a.a <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVa <- Ji
          
          sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) a = omegasq, 2) b = sigmasq
        
        sum.a.b <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVa <- Ji
          derVb <- Ri
          
          sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) a = omegasq, 2) c = tausq
        
        sum.a.c <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVa <- Ji
          #derVc <- Ii
          #sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))
          
          sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) b = sigmasq, 2) b = sigmasq
        
        sum.b.b <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVb <- Ri
          
          sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) b = sigmasq, 2) c = tausq
        
        sum.b.c <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVb <- Ri
          #derVc <- Ii
          #sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))
          
          sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) c = tausq, 2) c = tausq
        
        sum.c.c <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          #derVc <- Ii
          #sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))
          sum.c.c  <- sum.c.c + sum(Vi.inv * t(Vi.inv))
          
        }#for (i in 1 : nsubj)
        
        expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c,
                                            sum.a.b, sum.b.b, sum.b.c,
                                            sum.a.c, sum.b.c, sum.c.c), 
                                          ncol = 3, byrow = T)
        
        theta.new <- abs(as.numeric(theta.old - ginv(expected.hessian) %*% score.theta))
        
        if(silent == FALSE){
          
          cat("iteration = ", iter, "\n")
          cat("theta.old = ", theta.old, "\n")
          cat("beta.hat  = ", b.hat, "\n")
          cat("theta.new = ", theta.new, "\n")
          cat("sqrt.diff=", sqrt((theta.old - theta.new) %*% (theta.old - theta.new)), "\n")
          cat("score=", score.theta, "\n")      
          print("-----------------------------")
          
        }
        
        iter <- iter + 1
        
      }#while
      
      
      theta.old <- theta.new
      
      a <- theta.old[1]  ## omegasq
      b <- theta.old[2]  ## sigmasq
      c <- theta.old[3]  ## tausq
      
      ### betahat
      
      sum.left.beta <- sum.right.beta <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        
        Xi.transp <- t(Xi)
        
        sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
        sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi
        
      }#for (i in 1 : nsubj)
      
      b.hat    <- new_solve(sum.left.beta) %*% sum.right.beta
      b.varcov <- new_solve(sum.left.beta)
      
      ####### score for theta
      
      ### a = omegasq
      
      sum.a <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        ri       <- Yi - Xi %*% b.hat 
        
        derVa <- Ji
        
        sum.a <- sum.a + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVa))
        
      }#for (i in 1 : nsubj)
      
      
      ### b = sigmasq
      
      sum.b <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        ri       <- Yi - Xi %*% b.hat 
        
        derVb <- Ri
        
        sum.b <- sum.b + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVb))
        
      }#for (i in 1 : nsubj)
      
      
      # c = tausq
      
      sum.c <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        ri       <- Yi - Xi %*% b.hat 
        
        derVc <- Ii
        
        sum.c  <- sum.c + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVc))
        
      }#for (i in 1 : nsubj)
      
      score.theta <- 0.5 * matrix(c(sum.a, sum.b, sum.c))
      
      ######### INFORMATION MATRIX
      
      #### 1) a = omegasq, 2) a = omegasq
      
      sum.a.a <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        
        derVa <- Ji
        
        sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))
        
      }#for (i in 1 : nsubj)
      
      
      #### 1) a = omegasq, 2) b = sigmasq
      
      sum.a.b <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        
        derVa <- Ji
        derVb <- Ri
        
        sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))
        
      }#for (i in 1 : nsubj)
      
      
      #### 1) a = omegasq, 2) c = tausq
      
      sum.a.c <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        
        derVa <- Ji
        #derVc <- Ii
        #sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))
        
        sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv))
        
      }#for (i in 1 : nsubj)
      
      
      #### 1) b = sigmasq, 2) b = sigmasq
      
      sum.b.b <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        
        derVb <- Ri
        
        sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))
        
      }#for (i in 1 : nsubj)
      
      
      #### 1) b = sigmasq, 2) c = tausq
      
      sum.b.c <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        
        derVb <- Ri
        #derVc <- Ii
        #sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))
        
        sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv))
        
      }#for (i in 1 : nsubj)
      
      
      #### 1) c = tausq, 2) c = tausq
      
      sum.c.c <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        
        #derVc <- Ii
        #sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))
        
        sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% Vi.inv))
        
      }#for (i in 1 : nsubj)
      
      expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c,
                                          sum.a.b, sum.b.b, sum.b.c,
                                          sum.a.c, sum.b.c, sum.c.c), 
                                        ncol = 3, byrow = T)
      
      sd.theta <- sqrt(diag(ginv(-expected.hessian)))
      
      ## loglik
      
      sum.loglik <- 0
      
      for (i in 1 : nsubj){
        
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Xi    <- as.matrix(DM[[i]], nrow = ni)
        Yi    <- as.matrix(YM[[i]], nrow = ni)
        
        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
        Ii <- diag(ni)
        
        Vi <- a * Ji + b * Ri + c * Ii
        
        Vi.inv <- new_solve(Vi)
        
        ri     <- Yi - Xi %*% b.hat
        
        sum.loglik  <- sum.loglik + as.numeric(determinant(Vi, logarithm=TRUE)$modulus) +  t(ri) %*% Vi.inv %*% ri
        
      }#for (i in 1 : nsubj)
      
      max.loglik <- as.numeric(- 0.5 * ntotal * log(2 * pi) - 0.5 * sum.loglik)
      
      
      result <- cbind(c(c(b.hat), theta.new) , c(sqrt(diag(b.varcov)), sd.theta),
                      c(c(b.hat)/sqrt(diag(b.varcov)), NA, NA, NA))
      result <- cbind(result, 2 * (1 - pnorm(abs(result[, 3]))))
      colnames(result) <- c("Estimate", "Standard error", "Z-estimate", "p-value")
      rownames(result)[(nrow(result) - 2) : nrow(result)] <- c("omegasq", "sigmasq", "tausq")
      
      score.theta <- matrix(score.theta, nrow = 1)
      colnames(score.theta) <- c("omegasq", "sigmasq", "tausq")
      
      rand.varcov <- ginv(-expected.hessian) 
      colnames(rand.varcov) <- rownames(rand.varcov) <- c("omegasq", "sigmasq", "tausq")
      
      output           <- list()
      output$title     <- "Mixed effects model with random intercept and integrated Brownian motion"
      output$date      <- date()
      output$estimates <- result
      output$maxloglik <- max.loglik 
      output$score     <- score.theta
      output$fix.varcov  <- b.varcov 
      output$rand.varcov <- rand.varcov
      output
      
    }#ibm
    
    output
    
  }
