##
## Tests based on chisquare goodness-of-fit tests
##
## --------------------------------------------------------------------------

rvgt.chisq <- function (ftable)

  ## ------------------------------------------------------------------------
  ## Perform chisquare test on frequency table.
  ## ------------------------------------------------------------------------
  ## ftable : Object of class "rvgt.ftable" containing frequencies
  ## ------------------------------------------------------------------------
  ## Return:
  ## list of p values of cumulative frequencies
  ## ------------------------------------------------------------------------
{ 
  ## check arguments
  if (missing(ftable) || class(ftable) != "rvgt.ftable")
    stop ("Argument 'ftable' missing or invalid.")

  ## get table
  table <- ftable$count
  
  ## samplesize for one repetition
  n <- ftable$n

  ## number of repetitions
  r <- ftable$rep
  
  ## number of bins
  nbins <- ncol(table)

  ## probabilities under null hypothesis
  ubreaks <- ftable$ubreaks
  p0 <- diff(ubreaks)
  
  ## Vector to store p-values.
  pval <- numeric(r)
  
  ## array for computing cumulative frequencies
  fcum <- numeric(nbins)

  ## compute p-values of cumulative frequencies
  for (i in 1:r) {
    fcum <- fcum + table[i,]
    pval[i] <- chisq.test(fcum,p=p0)$p.value
  }
  
  ## return result as object of class "rvgt.htest"
  result <- list (type="chisq", n=n, rep=r, breaks=nbins+1, pval=pval)
  class(result) <- "rvgt.htest"

  return (result)    
}

## --------------------------------------------------------------------------

print.rvgt.htest <- function (x, ...) {
  cat("\nrvgtest - test:\n")
  cat("   type           =",x$type,"\n");
  cat("   sample size    =",x$n*x$rep,"\n");
  cat("   # break points =",x$breaks,"\n")
  cat("   p-value        =",x$pval[length(x$pval)],"\n\n")
}

## --------------------------------------------------------------------------
