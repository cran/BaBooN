
BBPMM.row <- function(misDataPat,
                      blockImp = length(misDataPat$blocks),
                      M=10,
                      outfile = NULL,
                      manWeights = NULL,
                      stepwise = TRUE,
                      verbose = TRUE,
                      tol=0.25,
                      setSeed = NULL,
                      ...)
{
  data.set <- misDataPat$data
  key <- misDataPat$key
  if(length(setdiff(blockImp, 1:length(misDataPat$blocks))) > 0) {
    stop(paste("blockImp =",as.character(setdiff(blockImp, 1:length(misDataPat$blocks))),
               "is not a subset of the number of",
         "different missing-data patterns (blocks)!\n"))}
  b.order <- order(blockImp)
  blockImp <- sort(blockImp)
  block <- misDataPat$blocks[blockImp]
  comp.names <- misDataPat$compNames
  n  <- nrow(data.set)
  R <- matrix(is.na(data.set), nrow=n)
  if (!is.null(key)) {
   pairlist <- vector(mode="list", length=M)
   names(pairlist) <- paste("M",1:M,sep="") 
  }
  weight.matrix <- vector(mode="list", length=M)
  names(weight.matrix) <- paste("M",1:M,sep="")
  model <- vector(mode="list", length=M)
  names(model) <- paste("M",1:M,sep="")
  y.hat <- vector(mode="list", length=M)
  names(y.hat) <- paste("M",1:M,sep="")
  impdata <- vector(mode="list", length=M)
  names(impdata) <- paste("M",1:M,sep="")
  miss <- function(x) {any(is.na(x)) }
  if (!is.null(key) && ncol(key) == 1) {
    Donid <- Recid <- names(key)[1]
  } else if (!is.null(key) && ncol(key) == 2) {
    Recid <- names(key)[1]
    Donid <- names(key)[2]
  } else if (is.null(key)) {
    pairlist <- NULL}
  varnames <- names(data.set)
  if (!is.null(setSeed)) set.seed  <- setSeed
  if(!is.null(manWeights)) {
    if(!is.vector(manWeights, mode = "list")) {
       manWeights <- list(manWeights)
       if (length(blockImp) > 1) {
       stop(paste("Only one vector with manual weights, but more than one",
                  "block specified for imputation!\n"))}
    } else if (is.vector(manWeights, mode = "list") & (length(blockImp) < length(manWeights))) {
      stop(paste("'manWeight' contains more elements than 'blockImp'!\n"))
    }
      if (any(unlist(manWeights) < 0)) {
        stop(paste("manWeights contains negative value(s)!\n"))}
    if (length(b.order) > 1) manWeights <- manWeights[b.order]
  }
  n <- nrow(data.set)
  l <- ncol(data.set)
  mis.pos <- vector(mode="list", length=length(block))
  obs.pos <- vector(mode="list", length=length(block))
  s.model <- vector(mode="list", length=length(block))
  nYC <- "YC"
  nM <- "M"
  while (any(varnames == nYC)) nYC <- paste(nYC,"C",sep="")
  while (any(varnames == nM)) nM <- paste(nM,"M",sep="")
  #########################################################################
  for (j in seq(along=block)) {
    mrow <- apply(as.matrix(data.set[ ,block[[j]]]), 1, miss)
    mis.pos[[j]] <- (1:n) [mrow == TRUE]
    obs.pos[[j]] <- (1:n) [mrow == FALSE]
    ## Test for available degress of freedom in the model
    mc.test <- qr(as.matrix(data.set[obs.pos[[j]], comp.names]),...)
    if((((length(obs.pos[[j]])-1) <= length(comp.names)) & stepwise == FALSE) |
       (mc.test$rank != length(comp.names) & stepwise == FALSE)) {
      stop("Block ",j, " has insufficient rank for imputation!\n")
    }
  }
  xvars <- paste(comp.names,collapse= ' + ') 
  ##########################################################################
  ### first loop for MI----------------------------------------------------
  for (m in 1:M) {
    if(!is.null(key)) {
      pairlist[[m]] <- vector(mode = "list", length=length(block))
      names(pairlist[[m]]) <- paste("block",seq(along=block),sep="")
    }
    weight.matrix[[m]] <- vector(mode = "list", length=length(block))
    names(weight.matrix[[m]]) <- paste("block",seq(along=block),sep="")
    model[[m]] <- vector(mode = "list", length=length(block))
    names(model[[m]]) <- paste("block",seq(along=block),sep="")
    y.hat[[m]] <- vector(mode = "list", length=length(block))
    names(y.hat[[m]]) <- paste("block",seq(along=block),sep="")
    impdata[[m]] <- data.set
    for (j in seq(along=block)) { # second loop for different blocks
      model[[m]][[j]] <- vector(mode="list", length=length(block[[j]]))
      names(model[[m]][[j]]) <- varnames[block[[j]]]
      S.xy <- NULL
      y.hat[[m]][[j]] <- matrix(nrow=n,ncol=length(block[[j]]))
      colnames(y.hat[[m]][[j]]) <- varnames[block[[j]]]
      co2 <- 0
      if(M > 1){
        BB.ind <- BayesBoot(ind.obs = obs.pos[[j]])
        BB.data <- data.set[BB.ind, ]
      }
      for (k in block[[j]]) {
        co2 <- co2+1
        s.model <- as.formula(ifelse(co2 == 1,
                                     paste(varnames[k],'~',xvars),
                                     paste(nYC,'~',xvars)))
        if (co2 == 1) {
          if(M == 1) {
            regmod <- lm(s.model, data=data.set, subset=obs.pos[[j]])
          } else if(M > 1) {
            BB.stab <- BB.mod.stab.glm(data=data.set,BB.data=BB.data,
                                       s.model=s.model)
            regmod <- BB.stab$model
            if (any(BB.stab$mislevpos == TRUE)) {
              warning(paste("Imputation ",m,", block ",j,
                            ": Bayesian Bootstrap dropped ",
                            "at least one category of a factor variable!\n",
                            sep=""))
            }
          }
          if (stepwise) regmod <- stepAIC(regmod, trace=0, k=log(n),
																						direction = "backward")
          var.T <- var(data.set[ ,k], na.rm=TRUE)
          var.U <- var(regmod$residuals)*length(obs.pos[[j]])/
          (length(obs.pos[[j]])-regmod$rank+1)
        } else if (co2 > 1) {
          xvars.e <- paste(varnames[block[[j]]][1:(co2-1)], collapse='+')
          s.model.e <- as.formula(paste(varnames[k],' ~ ',xvars.e))
          if(M == 1) {
            regmod.e <- lm(s.model.e, data=data.set, subset=obs.pos[[j]])
          } else if(M > 1) {
            regmod.e <- lm(s.model.e, data=BB.data)
          }
          YC <- numeric(n)
          YC[obs.pos[[j]]] <- regmod.e$residuals
          data.set <- as.data.frame(cbind(YC,data.set))
          names(data.set)[1] <- nYC
           if(M == 1) { 
            regmod <- lm(s.model, data=data.set, subset=obs.pos[[j]])
          } else if(M > 1) {
            BB.data <- as.data.frame(cbind(YC[BB.ind],BB.data))
            names(BB.data)[1] <- nYC
            BB.stab <- BB.mod.stab.glm(data=data.set,BB.data=BB.data,
                                       s.model=s.model)
            regmod <- BB.stab$model
          }
          if (stepwise) regmod <- stepAIC(regmod, trace=0, k=log(n),
																						direction = "backward")
          var.T <- var(YC, na.rm=TRUE)
          var.U <- var(regmod$residuals)*length(obs.pos[[j]])/
          (length(obs.pos[[j]])-regmod$rank+1)
        }
        model[[m]][[j]][[co2]] <- regmod 
        y.hat[[m]][[j]][ ,co2] <- predict(regmod,newdata=data.set)
        ## multicollinearity among ys
        if (is.na(var.T) || var.T < 1e-16) {
          S.xy[co2] <- 1e16
        } else {
          S.xy[co2] <- var.U}
        if (co2 > 1) { data.set <- data.set[ ,-1]
          if (M > 1) { BB.data <- BB.data[ ,-1] }}
        
      } ## end of loop k (incomplete variables)
      suppressWarnings(rm("BB.data"))
      gc()
      if (length(S.xy) > 1) {
        weight.matrix[[m]][[j]] <- diag(S.xy)
      } else {
        weight.matrix[[m]][[j]] <- S.xy}
      if (!is.null(manWeights) && length(manWeights[[j]]) > 0) {
        if (length(manWeights[[j]]) > 1) {
          weight.matrix[[m]][[j]] <- diag(manWeights[[j]]^(-1)*
          diag(weight.matrix[[m]][[j]]))
        } else if (length(manWeights[[j]]) == 1) {
          weight.matrix[[m]][[j]] <- weight.matrix[[m]][[j]]/manWeights[[j]]
        }
      }
      y.hat.obs <- y.hat[[m]][[j]][obs.pos[[j]], ]
      if (verbose) {
        cat(paste("Imputation ",m,": reciprocal weight matrix for block ",j,
                  ":\n",sep = ""))
        print(weight.matrix[[m]][[j]])
      }
      if (!is.null(key)) {
        pairlist[[m]][[j]] <- matrix(nrow=length(mis.pos[[j]]),ncol=2)}
      co3 <- 0
      for (i in mis.pos[[j]]) # third loop b) for the unobserved ys
        {
          co3 <- co3+1
          index <- obs.pos[[j]][apply(t(y.hat.obs),2,
                                 FUN = function(x) {
                                   t(y.hat[[m]][[j]][i, ] - x) %*%
                                     weight.matrix[[m]][[j]] %*%
                                       (y.hat[[m]][[j]][i, ] - x)})
                           == min(apply(t(y.hat.obs),2,
                                  FUN = function(x) {
                                    t(y.hat[[m]][[j]][i, ] - x) %*%
                                      weight.matrix[[m]][[j]] %*%
                                        (y.hat[[m]][[j]][i, ] - x)}))]
          if (length(index) > 1) {
            index <- sample(index, 1)
          } # random selection in case of several nearest neighbours
          if (!is.null(key)) {
            pairlist[[m]][[j]][co3, ] <- c(key[i,Recid],key[index,Donid]) }
          impdata[[m]][i, block[[j]]] <- data.set[index, block[[j]]]
        } ## end of i loop (missing values)
    } ## end of j loop (blocks)
    if (!is.null(key)) impdata[[m]] <- cbind(key, impdata[[m]])
  } ## end of m loop (MI)
  if (M == 1) impdata <- impdata[[m]]
  if (!is.null(outfile)) {
    outDAT <- matrix(nrow=n*M, ncol=l+1)
    if (M == 1) {
      outDAT <- cbind(impdata,1)
    } else {
      for (i in seq(along = impdata)) {
        outDAT[((i-1)*n+1):(i*n), ] <- cbind(as.matrix(impdata[[i]]),i)
      }
    }
    outDAT <- as.data.frame(outDAT)
    names(outDAT) <- c(varnames,nM) ## check for file ending
    if (grep(".",outfile) > 0) {
      lastDot <- max(which(strsplit(outfile,"")[[1]]=="."))
      StrL <- length(strsplit(outfile,"")[[1]])
      if ((StrL - lastDot) > 3) {
        outfile <- paste(outfile,".dat",sep="")
      }
    } else {
      outfile <- paste(outfile,".dat",sep="")
    }
    write.table(outDAT, file = outfile, sep = "\t",
                row.names = FALSE, quote = FALSE)
  }
  x <- list("impdata" = impdata, "weightMatrix" = weight.matrix, "model" = model,
            "pairlist" = pairlist, "indMatrix" = R)
  class(x) <- "imp"
  return(x)
}
