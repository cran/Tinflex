#############################################################################
##
##  Setup for generator.
##  (C only version)
##
#############################################################################

## See setup.C for R structures 'iv.struct' (for intervals) and
## 'gen.struct' for Tinflex object.

## --------------------------------------------------------------------------

Tinflex.setup.C <- function(lpdf, dlpdf, d2lpdf=NULL, ib, cT=0, rho=1.1, max.intervals=1001)
  ## -----------------------------------------------------------------------
  ## Setup: compute hat and squeeze for density.
  ## -----------------------------------------------------------------------
  ##   lpdf   ... log-density 
  ##   dlpdf  ... derivative of log-density
  ##   d2lpdf ... 2nd derivative of log-density
  ##   ib     ... interval boundaries of decomposition
  ##   cT     ... parameter for transformation (global or for each interval)
  ##   rho    ... performance parameter: requested upper bound for ratio
  ##              between area below hat and area below squeeze
  ##   max.intervals ... maximal numbers of intervals
  ## -----------------------------------------------------------------------
  ## return: S3 object of class "TinflexC" that contains
  ##   ivs      ... parameters for hat and squeeze for each interval
  ##                the last interval just stores the right boundary of domain
  ##   lpdf     ... log-density 
  ##   A.ht.tot ... total area below hat
  ##   A.sq.tot ... total area below squeeze
  ##   env      ... environment for evaluating log-density in C
  ##   iniv     ... string that describes initial intervals
  ##   Acum     ... cumulated areas  (=NULL; only for compatibilty with "Tinflex"
  ##   gt       ... guide table      (=NULL; only for compatibilty with "Tinflex"
  ##   Cgen     ... pointer to external C structure
  ## -----------------------------------------------------------------------
{
    ## Check arguments.
    if (!isTRUE(rho >= 1.0001)) {
        warning ("argument 'rho' too small or invalid, using default")
        rho <- 1.1
    }
    
    if (!isTRUE(max.intervals > 51)) {
        warning ("argument 'max.intervals' too small or invalid, using default")
        max.intervals <- 1001
    }
    
    if (missing(lpdf) || !is.function(lpdf)) {
        stop ("parameter 'lpdf' missing or invalid")
    }
    if (missing(dlpdf) || !is.function(dlpdf)) {
        stop ("parameter 'dlpdf' missing or invalid")
    }
    if (! (is.function(d2lpdf) || is.null(d2lpdf) ) ) {
        stop ("parameter 'd2lpdf' is invalid")
    }
    
    if (missing(ib)) {
        stop ("argument 'ib' is missing, with no default")
    }
    if (length(ib)<2 || length(ib) > max.intervals / 2) {
        stop ("argument 'ib' invalid")
    }
    
    ## The boundaries must be sorted.
    ib <- sort(ib)
    
    ## Check parameters for transformation.
    if (! ( is.numeric(cT) && (length(cT)==1 || length(cT)==length(ib)-1 )) ) {
        stop ("argument 'cT' invalid: its length must equal either 1 or number of intervals")
    }
    if ((! is.finite(ib[1]) && ! isTRUE(cT[1] > -1.)) ||
        (! is.finite(ib[length(ib)]) && ! isTRUE(cT[length(cT)] > -1.)) ) {
        stop ("(first and last) entry of argument 'cT' must be greater than -1 for unbounded domains")
    }
    
    ## ........................................................................
    
    ## Store some of the arguments for print.Tinflex().
    iniv <- paste(paste("    initial intervals ="),
                  paste(ib,collapse=" | "),
                  paste("\n                   cT ="),
                  paste(cT,collapse=" | "),
                  paste("\n                  rho =",rho,"\n"))

    ## We need an evironment for evaluating R expressions
    lpdf.env=parent.frame()

    ## Create S3 class that contains generator.
    generator <- list(
        ivs=NULL,          ## data for hat and squeeze
        lpdf=lpdf,         ## log-density
        A.ht.tot=NULL,     ## total area below hat
        A.sq.tot=NULL,     ## total area below hat
        env=lpdf.env,      ## environment for evaluating log-density
        iniv=iniv,         ## initial intervals (for print.Tinflex)
        Acum=NULL,         ## cumulated areas
        gt=NULL,           ## guide table
        haved2Tf= !is.null(d2lpdf), ## whether 2nd derivative is used (for print.Tinflex)
        Cgen=NULL          ## pointer to external C structure
    )
    class(generator) <- "TinflexC"

    ## Create TinflexC object
    generator$Cgen <- .Call(C_Tinflex_C_setup, generator,
                            lpdf.env, lpdf, dlpdf, d2lpdf,
                            ib, cT, rho, max.intervals);

    ## Return generator object. 
    return (generator)
} ## -- end of Tinflex.setup.C() -- ##

## --------------------------------------------------------------------------

Tinflex.sample.C <- function(gen, n=1) {
  ## ------------------------------------------------------------------------
  ## Draw a sample of size 'n'.
  ## (C version)
  ## ------------------------------------------------------------------------
  ##   gen ... S3 object generated by function 'Tinflex.setup'
  ##   n   ... sample size
  ## ------------------------------------------------------------------------
  ## Return: random sample.
  ## ------------------------------------------------------------------------

  if (class(gen)[1] != "TinflexC") {
      stop("Argument 'gen' is not of class 'TinflexC'")
  }

  .Call(C_Tinflex_C_sample, gen$Cgen, n)
}

## --------------------------------------------------------------------------

Tinflex.C2R <- function(genC) {
    ## ----------------------------------------------------------------------
    ## Convert object of class 'TinflexC' into 'Tinflex'
    ## ----------------------------------------------------------------------
    ##   gen ... S3 object generated by function 'Tinflex.setup.C'
    ## ----------------------------------------------------------------------
    ## Return: object of class 'Tinflex'
    ## ----------------------------------------------------------------------

    if (class(genC)[1] != "TinflexC") {
        stop("Argument 'gen' is not of class 'TinflexC'")
    }

    ## get list with data from 'genC'
    Cdata <- .Call(C_Tinflex_C_2_R, genC$Cgen);

    ## copy all data into new object
    gen <- genC
    gen$A.ht.tot <- Cdata$A.ht.tot
    gen$A.sq.tot <- Cdata$A.sq.tot
    gen$Acum <- Cdata$Acum
    gen$gt <- Cdata$gt

    ## data for intervals
    ivs <- Cdata$ivs
    rownames(ivs) <- iv.struct
    ## the following data are irrelevant but changed to NA
    ## as this is the entry in an object of class 'Tinflex'
    NA.fields <- c("ht.a","ht.b","ht.y","sq.a","sq.b","sq.y","A.ht","A.sq","type")
    ivs[NA.fields,dim(ivs)[2]] <- NA
    ivs["next",dim(ivs)[2]] <- Inf
    gen$ivs <- ivs

    ## remove obsolete pointer
    gen$Cgen <- NULL
    
    ## set class name
    class(gen) <- "Tinflex"

    ## return
    gen
    
} ## -- end of Tinflex.C2R() -- ##

## --------------------------------------------------------------------------
