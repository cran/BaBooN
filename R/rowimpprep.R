                                   


#################################################################
### preparation for row-wise imputation

rowimpPrep <- function(data,
                       ID=NULL,
                       verbose=TRUE)
{
  data <- as.data.frame(data)
  key <- NULL
  if (length(ID) > 2) stop("argument 'ID' contains more than two elements!")
  if (!is.null(ID)) {
    if (is.character(ID)) ID <- which(is.element(names(data),ID))
    key <- as.data.frame(data[ ,ID])
    data <- data[ ,-ID] } else {key <- NULL}
  block <- list()
  ##org.names <- var.names <- names(data) # variable names
  l <- ncol(data) # no. variables
  miss <- function(x) {any(is.na(x)) }
  mvar <- apply(data, 2, miss)  #which columns have missings values?
  p.impvar <- p.impvar2 <- (1:l)[mvar] #position of columns with missing
  ## values
  ## constant regressors?
  comp <- data[ ,-p.impvar] # completely observed variables 
  ## values
  comp.names <- names(comp) # names of completely observed variables
  dat.isna <- is.na(data)  # indicator matrix for missing values
  co1 <- 0
  for (i in 1:length(p.impvar))
    ## creates list with positions for identical columns (blocks)
    {
      co1 <- co1+1
      if(length(p.impvar2) > 1)
        {
          aa <- (dat.isna[,p.impvar2[1]] == dat.isna[,p.impvar2]) - 1
          identical <- which(apply(aa, 2, sum) == 0)
          ## which columns are identical to the first column with missing
          ## values
          block[[co1]] <- p.impvar2[identical]
        } else if (length(p.impvar2) == 1){
          block[[co1]] <- p.impvar2
          break  # breaks off if one column is left
        } else {
          break  # breaks off if no column is left
        }
      p.impvar2 <- p.impvar2[-identical]
    }
  block.names <- lapply(block,FUN=function (x) names(data[x]))
  names(block.names) <- names(block) <- paste("block",1:length(block))
  if (verbose) {
    cat("variables with missing values ordered by identical patterns:\n")
    print(block.names)
  }
  x <- list("data"=data, "key"=key, "blocks"=block, "blockNames"=block.names,
              "compNames"=comp.names)
  class(x) <- "impprep"
  return(x)
}
