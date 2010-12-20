

BBPMM <- function(data,
                  M = 10,
                  nIter = 10,
                  outfile = NULL,
                  ignore = NULL,
                  vartype = NULL,
                  stepwise = TRUE,
                  maxit = 3,
                  verbose = TRUE,
                  setSeed = NULL,
                  ...)
  {
    DAT <- as.data.frame(data)
    orgnames <- varnames <- names(DAT)
		n <- nrow(DAT) 
    l <- ncol(DAT)
    if (!is.null(ignore)) {
      if(length(setdiff(ignore,varnames) > 0) | length(setdiff(ignore,1:l) > 0)) {
        stop("'ignore' is no subset of either number of columns or variable names!\n")
      }
      if (is.character(ignore)) {
        ig.pos <- is.element(varnames, ignore)
        ignore <- which(ig.pos)
      } else {
        ig.pos <- is.element(1:l, ignore)
      }
      not.inc <- as.data.frame(DAT[ ,ignore])
      varnames <- varnames[-ignore] 
      DAT <- DAT[ ,-ignore]
      l <- ncol(DAT)
    }
    if (l <= 2) {
      ## add up to two columns with random noise to ensure variability for PMM
      DAT <- as.data.frame(cbind(DAT,matrix(runif((3-l)*n),nrow=n)))
      varnames <- names(DAT)
      onevar <- TRUE
    }
    ## take over class from data.frame or administer classes
    if (!is.null(vartype)) {
      if (length(vartype) != l) {
        stop("Error: Number of flagged variables in 'vartype'",
             "does not match number of (remaining) variables in",
             "data set!\n")
      } else if (any(vartype != "C" & vartype != "M")) {
        stop("Error: 'vartype' contains wrong character(s)!\n")
      }
      f.pos <- which(vartype == "C")
      if (length(f.pos) > 0) {
        DAT[,f.pos] <- lapply(DAT[,f.pos], as.factor)}
      m.pos <- which(vartype == "M")
      if (length(m.pos) > 0) {
        DAT[,m.pos] <- lapply(DAT[,m.pos], as.numeric)}
    }
    M.DAT <- vector(M, mode = "list")
    if (!is.null(setSeed)) set.seed  <- setSeed
    nM <- "M"
    while (any(varnames == nM)) nM <- paste(nM,"M",sep="")
    onevar <- FALSE
    org.l <- ncol(DAT)
    ## indicator matrix for missing values
    R <- matrix(is.na(DAT), nrow=n)
    mis.num <- colSums(is.na(DAT))
    mis.overview <- paste("number of missing values ", names(DAT),": ",
                          mis.num, sep="")
    if (verbose) print(mis.overview)
    ## new variable order
    n.order <- order(mis.num)
    o.order <- order(n.order)
    DAT <- as.data.frame(DAT[ ,n.order])
    varnames <- varnames[n.order]
    mvar <- apply(DAT, 2, FUN = function(x) {any(is.na(x)) })
    p.impvar <- (1:l)[mvar]
    p.comp <- (1:l)[!mvar]
    i.mis <- vector(length(p.impvar), mode = "list")
    i.obs <- vector(length(p.impvar), mode = "list")
    for (j in seq(along = p.impvar)) {
      i.mis[[j]] <- (1:n)[is.na(DAT[ ,p.impvar[j]])]
      i.obs[[j]] <- (1:n)[!is.na(DAT[ ,p.impvar[j]])]
    }
    ## starting solution
    startSol <- TRUE
    ##+++++++++++++++++++++++++PMM++++++++++++++++++++++++++++++++++++
    ## Sequential Regression with Predictive Mean Matching
    for (m in 1:M) {
      co <- 0
      iterate <- TRUE
      while (iterate) {
        ##--------------first loop for iterations-----------------------
        if (!startSol) co <- co + 1
        co2 <- 0
        if (verbose & !startSol) {
          cat(paste("Imputation ", m," of ",M ,": iteration ", co,
                    sep=""), "\n") }
        ##-------------- Bayesian Bootstrap --------------------------
        if (M > 1) {
          ind1 <- BayesBoot(ind.obs = 1:n)
          ## Bayesian Bootstrap: draw n times with replacement as basis for
          ## imputation model parameter estimates
        }
        if (verbose & !startSol) cat("Variable:")
        ##----second loop for every variable with missing values------
        for (j in p.impvar) {
          if (verbose & !startSol) cat("",varnames[j])
          if (verbose & !startSol & j == rev(p.impvar)[1]) cat("\n ")
          co2 <- co2 + 1
          if (startSol & length(p.comp) == 0) {
            ## Hotdeck imputation for starting solution if no variable
            ## completely observed
            DAT[i.mis[[co2]],j] <- sample(DAT[i.obs[[co2]],j],
                                             length(i.mis[[co2]]),
																						 replace = TRUE)
            p.comp <- j
          }
          y <- DAT[ ,j]
          if(startSol){
            xvars <- paste(c(varnames[p.comp],varnames[p.impvar[0:(co2-1)]]),
                           collapse=' + ')
          } else {
            xvars <- paste(varnames[-j], collapse = ' + ')
          }
          s.model <- as.formula(paste(varnames[j],'~',xvars))
          if (is.numeric(y)) {
            if (M == 1 | startSol) {
              regmod <- lm(s.model, data=DAT, subset = i.obs[[co2]])
            } else if (M > 1 & !startSol) {
							BB.data <- DAT[ind1, ]
              BB.stab <- BB.mod.stab.glm(data=DAT, BB.data=BB.data,
                                         s.model=s.model)
              regmod <- BB.stab$model}
            if (stepwise) regmod <- stepAIC(regmod, trace=0, k=log(n),
																						direction = "backward")
            y.pred <- predict(regmod, newdata=DAT)
            y.pred.mis <- y.pred[i.mis[[co2]]]
            y.pred.obs <- y.pred[i.obs[[co2]]]
            nextlist <- numeric(length(y.pred.mis))
            for (i in seq(along=y.pred.mis)) {
              nextlist[i] <- PMMsearchMet(yHatMis = y.pred.mis[i],
                                          yHatObs = y.pred.obs)
            }
          } else if (is.factor(y) & length(table(y)) > 2) {
            if (M == 1 | startSol) {
              options(warn = -1)
              regmod <- multinom(s.model,data=DAT,trace=F,
                                 subset = i.obs[[co2]])
              options(warn = 0)
            } else if ( M > 1 & !startSol) {
							BB.data <- DAT[ind1, ]
              BB.stab <- BB.mod.stab.mlog(data=DAT,
                                          BB.data=BB.data,
                                          s.model=s.model)
              regmod <- BB.stab$model
            }
            if (stepwise) regmod <- stepAIC(regmod, trace=0, k=log(n),
																						direction = "backward")
            y.pred <- predict(regmod, newdata=DAT, type="probs")
            y.pred[y.pred > 0.999] <- 0.999
            y.pred[y.pred < 0.001] <- 0.001
            l.y.pred <- log(y.pred/(1-y.pred))
            y.pred.mis <- l.y.pred[i.mis[[co2]], ]
            y.pred.obs <- l.y.pred[i.obs[[co2]], ]
            ## calculate outer product for all obs/mis columns 
            m.dist <- matrix(rep(0,nrow(y.pred.mis)*nrow(y.pred.obs)),
                             nrow=nrow(y.pred.mis))
            for (j in 1:ncol(y.pred)){
              m.dist <- m.dist+outer(y.pred.mis[ ,j],y.pred.obs[ ,j],FUN="-")^2
            }
            nextlist <- max.col(as.matrix(m.dist*(-1)),
                                ties.method="random")
            rm("m.dist")
          } else if (is.factor(y) & length(table(y)) == 2) {
            if (M == 1 | startSol) {
              regmod <- glm(s.model, data=DAT,
                            family = binomial(link="logit"),
                            subset = i.obs[[co2]])
            } else if (M > 1 & !startSol) {
							BB.data <- DAT[ind1, ]
              BB.stab <- BB.mod.stab.glm(data=DAT, BB.data=BB.data,
                                         s.model=s.model, model="binomial")
              regmod <- BB.stab$model
            }
            if (stepwise) regmod <- stepAIC(regmod, trace=0, k=log(n),
																						direction = "backward")
            y.pred <- predict(regmod, newdata=DAT)
            y.pred.mis <- y.pred[i.mis[[co2]]]
            y.pred.obs <- y.pred[i.obs[[co2]]]
            nextlist <- numeric(length(y.pred.mis))
            for (i in seq(along=y.pred.mis)) {
              nextlist[i] <- PMMsearchMet(yHatMis = y.pred.mis[i],
                                          yHatObs = y.pred.obs)
            }
          }
          DAT[i.mis[[co2]],j] <- y[i.obs[[co2]]][nextlist]
          gc()
        }
        startSol <- FALSE
        if (co == nIter) iterate <- FALSE
      } ## end of iterate cycle
      if (onevar == TRUE) {
        M.DAT[[m]] <- as.data.frame(DAT[,1:l])
        names(M.DAT[[m]]) <- orgnames
      }
      if (is.null(ignore)) {
        M.DAT[[m]] <- DAT[ ,o.order]
      } else {
        M.DAT[[m]] <- as.data.frame(matrix(nrow=n,ncol=org.l))
        M.DAT[[m]][ ,!ig.pos] <- DAT[ ,o.order]
        M.DAT[[m]][ ,ig.pos] <- not.inc
        names(M.DAT[[m]]) <- orgnames
      }
      if (M == 1) M.DAT <- M.DAT[[m]]
    } ## end of m cylce
    if (!is.null(outfile)) {
      outDAT <- matrix(nrow=n*M, ncol=l+1)
      if (M == 1) {
        outDAT <- cbind(M.DAT,1)
      } else {
        for (i in seq(along = M.DAT)) {
          outDAT[((i-1)*n+1):(i*n), ] <- cbind(as.matrix(M.DAT[[i]]),i)
        }
      }
      outDAT <- as.data.frame(outDAT)
      names(outDAT) <- c(orgnames,nM)
      ## check for file ending
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
    x <- list("impdata"=M.DAT, "misOverview"=mis.overview, "indMatrix"=R,
							"M"=M, "nIter"=nIter)
    class(x) <- "imp"
    return(x)
  }

