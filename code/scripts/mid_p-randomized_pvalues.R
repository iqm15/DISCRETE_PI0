# function for computing mid-p and randomized p-values 

#-------------------------------------------------------
#      p-value 'matching' function
#-------------------------------------------------------
match.pvals <- function(pCDFlist, raw.pvalues){

  m <- length(raw.pvalues)
  pvec <- raw.pvalues
  in.CDF <- numeric(m)
  for (k in (1:m)) {
    in.CDF[k] <- match(pvec[k], pCDFlist[[k]])
    if (is.na(in.CDF[k])){
      in.CDF[k] <- which.min(abs(pCDFlist[[k]] - pvec[k]))
      pvec[k] <- pCDFlist[[k]][in.CDF[k]]
      warning("Since the p-value ", raw.pvalues[k], " is not a value of the CDF of the p-value,\n                  the p-value is rounded to be ", 
              pCDFlist[[k]][in.CDF[k]], call. = F)
    }
  }
  return(pvec)
}



#------------------------------------------------------------------------------------------
#      Compute mid-p values
#------------------------------------------------------------------------------------------
midp.values <- function(pCDFlist, raw.pvalues){
 
  m <- as.integer(length(raw.pvalues))
  # stepf <- lapply(pCDFlist, function(x) stepfun(x, c(0, x)))
  pvec <- match.pvals(pCDFlist, raw.pvalues)
  
  midpCDFproba <- list()
  midpsupport <- list()
  
  obs_midp <- numeric(m)
  for(k in (1:m)){
    # compute the mid-pvalue
    idx <- match(pvec[k], unique(pCDFlist[[k]]))
    y <- c(0, unique(pCDFlist[[k]]))
    y.diff <- diff(y)
    obs_midp[k] <- pvec[k] - 0.5 * y.diff[idx]
    
    # get the "truncated" support and the new probabilities
    
    # midpsupport[[k]] <- pCDFlist[[k]][which((1 - duplicated(stepf[[k]](pCDFlist[[k]] + 0.5 * y.diff[idx]))) == 1)]
    # midpCDFproba[[k]] <- diff(c(0, unique(stepf[[k]](pCDFlist[[k]] + 0.5 * y.diff[idx]))))
    
    midpsupport[[k]] <- unique(pCDFlist[[k]]) - 0.5 * y.diff
   
  }
  
  output <- list(obs_midp = obs_midp, midpsupport = midpsupport)
  return(output)
}


#------------------------------------------------------------------------------------------
#      Compute randomized values
#------------------------------------------------------------------------------------------
randomizedp.values <- function(pCDFlist, raw.pvalues){
  
  m <- as.integer(length(raw.pvalues))
  
  pvec <- match.pvals(pCDFlist, raw.pvalues)
  
  mpv <- numeric(m)
  for(k in (1:m)){
    u <- runif(1)
    idx <- match(pvec[k], pCDFlist[[k]])
    y <- c(0, pCDFlist[[k]])
    y.diff <- diff(y)
    mpv[k] <- pvec[k] - (1 - u) * y.diff[idx]
  }
  return(mpv)
}
