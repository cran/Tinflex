#############################################################################
##
##  Setup for generator.
##
#############################################################################

## We store all data about the hat and squeeze in particular intervals
## in a matrix where each interval is stored in a column.
## We need the following parameter for each interval:

iv.struct <- 
    c(
        ## Data for hat and squeeze:
        "x",        ## left boundary of interval
        "c",        ## parameter for transformation
        "ht.a",     ## intercept for hat in transformed scale
        "ht.b",     ## slope for hat in transformed scale
        "ht.y",     ## anchor point of linear function
        "sq.a",     ## intercept for squeeze in transformed scale
        "sq.b",     ## slope for squeeze in transformed scale
        "sq.y",     ## anchor point of linear function
        "A.ht",     ## area below hat
        "A.sq",     ## area below squeeze
        ## Data required for setup only:
        "type",     ## type of interval
        "Tfx",      ## log-density at left boundary 
        "dTfx",     ## derivative of log-density at left boundary 
        "d2Tfx",    ## 2nd derivative of log-density at left boundary
        "next"      ## poor man's linked list:
        ## index of next interval to the right ('Inf' for last interval)
    )

## Names of list elements for generator: 

gen.struct <-
    c(
        "ivs",      ## data for hat and squeeze
        "lpdf",     ## log-density
        "A.ht.tot", ## total area below hat
        "A.sq.tot", ## total area below squeeze
        "env",      ## environment for evaluating log-density
        "iniv",     ## initial intervals (for print.Tinflex)
        "Acum",     ## cumulated areas
        "gt",       ## guide table
        "haved2Tf"  ## whether 2nd derivative is used (for print.Tinflex)
    )

## --------------------------------------------------------------------------

Tinflex.setup <- function(lpdf, dlpdf, d2lpdf=NULL, ib, cT=0, rho=1.1, max.intervals=1001) {
    ## ----------------------------------------------------------------------
    ## Setup: compute hat and squeeze for density.
    ## ----------------------------------------------------------------------
    ##   lpdf   ... log-density 
    ##   dlpdf  ... derivative of log-density
    ##   d2lpdf ... 2nd derivative of log-density
    ##   ib     ... interval boundaries of decomposition
    ##   cT     ... parameter for transformation (global or for each interval)
    ##   rho    ... performance parameter: requested upper bound for ratio
    ##              between area below hat and area below squeeze
    ##   max.intervals ... maximal numbers of intervals
    ## ----------------------------------------------------------------------
    ## return: S3 object of class "Tinflex" that contains
    ##   ivs      ... parameters for hat and squeeze for each interval
    ##                the last interval just stores the right boundary of domain
    ##   lpdf     ... log-density 
    ##   A.ht.tot ... total area below hat
    ##   A.sq.tot ... total area below squeeze
    ##   env      ... environment for evaluating log-density in C
    ##   iniv     ... string that describes initial intervals
    ##   Acum     ... cumulated areas
    ##   gt       ... guide table
    ##   haved2Tf ... whether 2nd derivative is used (for print.Tinflex)
    ## ----------------------------------------------------------------------

    ## Check arguments.
    if (!isTRUE(rho >= 1.0001)) {
        warning ("argument 'rho' too small or invalid, using 1.0001")
        rho <- 1.0001
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

    if (is.null(d2lpdf)) {
        generator <- setup.no2ndD  (lpdf, dlpdf,         ib, cT, rho, max.intervals)
    }
    else {
        generator <- setup.with2ndD(lpdf, dlpdf, d2lpdf, ib, cT, rho, max.intervals)
    }

    ## return
    generator
    
} ## -- end of Tinflex.setup() -- ##

## --------------------------------------------------------------------------

setup.with2ndD <- function(lpdf, dlpdf, d2lpdf, ib, cT, rho, max.intervals) {
    ## ----------------------------------------------------------------------
    ## Setup: compute hat and squeeze for density.
    ## Version that makes use of 2nd derivative.
    ## ----------------------------------------------------------------------
    ##   lpdf   ... log-density 
    ##   dlpdf  ... derivative of log-density
    ##   d2lpdf ... 2nd derivative of log-density
    ##   ib     ... interval boundaries of decomposition
    ##   cT     ... parameter for transformation (global or for each interval)
    ##   rho    ... performance parameter: requested upper bound for ratio
    ##              between area below hat and area below squeeze
    ##   max.intervals ... maximal numbers of intervals
    ## ----------------------------------------------------------------------
    ## return: S3 object of class "Tinflex" that contains
    ##   ivs      ... parameters for hat and squeeze for each interval
    ##                the last interval just stores the right boundary of domain
    ##   lpdf     ... log-density 
    ##   A.ht.tot ... total area below hat
    ##   A.sq.tot ... total area below squeeze
    ##   env      ... environment for evaluating log-density in C
    ##   iniv     ... string that describes initial intervals
    ##   Acum     ... cumulated areas
    ##   gt       ... guide table
    ## ----------------------------------------------------------------------

    ## Store some of the arguments for print.Tinflex().
    iniv <- paste(paste("    initial intervals ="),
                  paste(ib,collapse=" | "),
                  paste("\n                   cT ="),
                  paste(cT,collapse=" | "),
                  paste("\n                  rho =",rho,"\n"))

    ## Get parameter 'cT' for transformation.
    if (length(cT)==1) {
        cT <- rep(cT[1], length(ib))
    }
    else { ## Just duplicate last element.
        cT <- append(cT, cT[length(cT)])
    }
    
    ## Create a table for hat and squeeze.
    params <- matrix(nrow=length(iv.struct), ncol=max.intervals)
    rownames(params) <- iv.struct
    
    ## We need transformed densities and their derivatives on
    ## both boundaries of each interval.
    ## However, the boundary to the right of one interval is the
    ## boundary to the left of the consecutive interval.
    ## If their parameters 'c' coincide we only need to store the
    ## transformed density for the boundary to the left.
    ## If these parameters differ we insert an interval of length 0
    ## in order to store the transformed density for the boundary
    ## to the right. 
    
    n.ivs <- 0
    for (i in 1:length(ib)) {
        n.ivs <- n.ivs+1
        params["x", n.ivs] <- ib[i]
        params["c", n.ivs] <- cT[i]
        params[c("Tfx","dTfx","d2Tfx"), n.ivs] <- Tfdd(lpdf, dlpdf, d2lpdf, cT[i], ib[i])
        params["next", n.ivs] <- n.ivs+1
        
        ## Parameter 'cT' changes for next interval.
        ## Thus we insert an interval of length 0.
        if (isTRUE(i < length(ib)) && !identical(cT[i], cT[i+1])) {
            n.ivs <- n.ivs+1
            params["x", n.ivs] <- ib[i+1]
            params["c", n.ivs] <- cT[i]
            params[c("Tfx","dTfx","d2Tfx"), n.ivs] <- Tfdd(lpdf, dlpdf, d2lpdf, cT[i], ib[i+1])
            params["next", n.ivs] <- n.ivs+1
        }
    }
    
    ## Update number of intervals (after creating initial intervals).
    n.ivs <- n.ivs-1
    
    ## Terminate linked list.
    params["next", n.ivs+1] <- Inf

    ## NOTE: params[,n.ivs+1] now stores just the right boundary of the domain. 

    ## Compute parameters for hat and squeeze for initial intervals.
    for (i in 1:n.ivs) {
        params["type",i] <- type.with2ndD(left=params[,i], right=params[,i+1])
        params[,i] <- hat.iv(left=params[,i], right=params[,i+1])
        ## DEBUG:
        ## print(params[,i])
    }
    ## DEBUG:
    ## print(params[,n.ivs+1])
    
    ## Compute total areas for initial hat and squeeze.
    A.ht.tot <- sum(params["A.ht", 1:n.ivs])
    A.sq.tot <- sum(params["A.sq", 1:n.ivs])
    
    ## We have to split intervals where the area
    ## between hat and squueze is too large.
    while (! isTRUE(A.ht.tot / A.sq.tot <= rho)) {
        
        ## Compute average area.
        threshold <- 0.99 * (A.ht.tot - A.sq.tot) / n.ivs
        
        ## Find all intervals where area between hat and squeeze is too large.
        toolarge <- which(params["A.ht", 1:(n.ivs+1)] - params["A.sq", 1:(n.ivs+1)] >= threshold)
        if (length(toolarge) <= 0) break;
        if (n.ivs + length(toolarge) >= max.intervals) {
            warning ("maximum number of intervals exceeded")
            break;
        }
        
        ## Split all intervals.
        for (i in toolarge) {
            ## Get boundary points of interval.
            left <- params[,i]
            right <- params[,params["next",i]]
            
            ## Splitting point for interval (use "arc mean").
            p <- arc.mean(left["x"],right["x"])
            
            ## Increment number of intervals
            n.ivs <- n.ivs+1
            
            ## Create new interval on r.h.s.
            params["x", n.ivs+1] <- p
            params["c", n.ivs+1] <- left["c"]
            params[c("Tfx","dTfx","d2Tfx"), n.ivs+1] <- Tfdd(lpdf, dlpdf, d2lpdf, left["c"], p)
            params["type",n.ivs+1] <- type.with2ndD(left=params[,n.ivs+1], right=right)
            params[,n.ivs+1] <- hat.iv(left=params[,n.ivs+1], right)
            params["next",n.ivs+1] <- left["next"]
            
            ## Update interval on l.h.s.
            params["type",i] <- type.with2ndD(left=params[,i], right=params[,n.ivs+1])
            params[,i] <- hat.iv(left=params[,i], params[,n.ivs+1])
            params["next",i] <- n.ivs+1
        }
        
        ## Update total areas.
        A.ht.tot <- sum(params["A.ht",], na.rm=TRUE)
        A.sq.tot <- sum(params["A.sq",], na.rm=TRUE)
    }
    
    ## Check result.
    if (! is.finite(A.ht.tot)) {
        stop ("Cannot create hat function. A_hat is not finite: A_hat=", A.ht.tot)
    } 
    
    if (! isTRUE(A.ht.tot > 0.)) {
        stop ("Cannot create hat function. A_hat is not > 0: A_hat=", A.ht.tot)
    } 
    
    if (! is.finite(A.sq.tot)) {
        stop ("Cannot create squeeze. A_squeeze is not finite: A_squeeze=", A.sq.tot)
    } 
    
    if (! isTRUE(A.sq.tot >= 0.)) {
        stop ("Cannot create squeeze. A_squeeze is not >= 0: A_squeeze=", A.sq.tot)
    } 
    
    if (isTRUE(A.ht.tot < A.sq.tot))
        stop ("Invalid input: A_hat < A_squeeze!")
    
    if (!isTRUE(A.ht.tot / A.sq.tot <= rho))
        warning ("ratio A_hat / A_squeeze = ", A.ht.tot / A.sq.tot,
                 " is larger than rho = ",rho)
    
    ## Truncate working array.
    params <- params[,1:(n.ivs+1)]
    
    ## Sort array (so that the linked list becomes obsolete).
    params <- params[,order(params["x",])]
    
    ## Compute guide table.
    ## For R version use
    ##   guidetable <- make.guidetable(params)
    ##   Acum <- guidetable$Acum
    ##   gt <- guidetable$gt
    ## C version:
    Acum <- numeric(n.ivs)
    gt <- integer(n.ivs)
    A.ht.tot <- .Call("Tinflex_RC_make_guide_table", params, Acum, gt)
    
    ## Create S3 class that contains generator.
    generator <- list(
        ivs=params,          ## data for hat and squeeze
        lpdf=lpdf,           ## log-density
        A.ht.tot=A.ht.tot,   ## total area below hat
        A.sq.tot=A.sq.tot,   ## total area below hat
        env=parent.frame(2), ## environment for evaluating log-density
        iniv=iniv,           ## initial intervals (for print.Tinflex)
        Acum=Acum,           ## cumulated areas
        gt=gt,               ## guide table
        haved2Tf=TRUE        ## whether 2nd derivative is used (for print.Tinflex)
    )
    class(generator) <- "Tinflex"
    
    ## Return generator object. 
    generator

} ## -- end of setup.with2ndD() -- ##

## --------------------------------------------------------------------------

setup.no2ndD <- function(lpdf, dlpdf, ib, cT=0, rho=1.1, max.intervals=1001) {
    ## ----------------------------------------------------------------------
    ## Setup: compute hat and squeeze for density.
    ## Version without 2nd derivative of transformed density.
    ## ----------------------------------------------------------------------
    ##   lpdf   ... log-density 
    ##   dlpdf  ... derivative of log-density
    ##   ib     ... interval boundaries of decomposition
    ##   cT     ... parameter for transformation (global or for each interval)
    ##   rho    ... performance parameter: requested upper bound for ratio
    ##              between area below hat and area below squeeze
    ##   max.intervals ... maximal numbers of intervals
    ## ----------------------------------------------------------------------
    ## return: S3 object of class "Tinflex" that contains
    ##   ivs      ... parameters for hat and squeeze for each interval
    ##                the last interval just stores the right boundary of domain
    ##   lpdf     ... log-density 
    ##   A.ht.tot ... total area below hat
    ##   A.sq.tot ... total area below squeeze
    ##   env      ... environment for evaluating log-density in C
    ##   iniv     ... string that describes initial intervals
    ##   Acum     ... cumulated areas
    ##   gt       ... guide table
    ## ----------------------------------------------------------------------

    ## ......................................................................
    
    ## Store some of the arguments for print.Tinflex().
    iniv <- paste(paste("    initial intervals ="),
                  paste(ib,collapse=" | "),
                  paste("\n                   cT ="),
                  paste(cT,collapse=" | "),
                  paste("\n                  rho =",rho,"\n"))
    
    ## Get parameter 'cT' for transformation.
    if (length(cT)==1) {
        cT <- rep(cT[1], length(ib))
    }
    else { ## Just duplicate last element.
        cT <- append(cT, cT[length(cT)])
    }
    
    ## Create a table for hat and squeeze.
    params <- matrix(nrow=length(iv.struct), ncol=max.intervals)
    rownames(params) <- iv.struct
    
    ## We need transformed densities and their derivatives on
    ## both boundaries of each interval.
    ## However, the boundary to the right of one interval is the
    ## boundary to the left of the consecutive interval.
    ## If their parameters 'c' coincide we only need to store the
    ## transformed density for the boundary to the left.
    ## If these parameters differ we insert an interval of length 0
    ## in order to store the transformed density for the boundary
    ## to the right. 

    n.ivs <- 0
    for (i in 1:length(ib)) {
        n.ivs <- n.ivs+1
        params["x", n.ivs] <- ib[i]
        params["c", n.ivs] <- cT[i]
        params[c("Tfx","dTfx","d2Tfx"), n.ivs] <- Tfdd(lpdf, dlpdf, NULL, cT[i], ib[i])
        params["next", n.ivs] <- n.ivs+1
        
        ## Parameter 'cT' changes for next interval.
        ## Thus we insert an interval of length 0.
        if (isTRUE(i < length(ib)) && !identical(cT[i], cT[i+1])) {
            n.ivs <- n.ivs+1
            params["x", n.ivs] <- ib[i+1]
            params["c", n.ivs] <- cT[i]
            params[c("Tfx","dTfx","d2Tfx"), n.ivs] <- Tfdd(lpdf, dlpdf, NULL, cT[i], ib[i+1])
            params["next", n.ivs] <- n.ivs+1
        }
    }

    ## Update number of intervals (after creating initial intervals).
    n.ivs <- n.ivs-1
    
    ## Terminate linked list.
    params["next", n.ivs+1] <- Inf

    ## NOTE: params[,n.ivs+1] now stores just the right boundary of the domain. 

    ## Determine type of initial intervals.
    for (i in 1:n.ivs) {
        params["type",i] <- type <-
            type.no2ndD.init(left=params[,i], right=params[,i+1], lpdf, dlpdf)
        
        if (type == IIb_IVa__IIa_IVa ||       ## Case (3.3.3)
            type == IIIa_IVb__IIIb_IVb  ) {   ## Case (4.3.3)
            ## we have to split the interval

            ## Splitting point for interval (use "arc mean").
            left <- params[,i]
            p <- arc.mean(left["x"],params["x",i+1])

            ## Increment number of intervals
            n.ivs <- n.ivs+1
            
            ## Create new interval for r.h.s.
            params["x", n.ivs+1] <- p
            params["c", n.ivs+1] <- left["c"]
            params[c("Tfx","dTfx","d2Tfx"), n.ivs+1] <- Tfdd(lpdf, dlpdf, NULL, left["c"], p)
            params["next",n.ivs+1] <- left["next"]
            params["type",n.ivs+1] <- if (type == IIb_IVa__IIa_IVa) { IIa_IVa } else { IIIb_IVb }
            
            ## Update interval on l.h.s.
            params["next",i] <- n.ivs+1
            params["type",i] <- if (type == IIb_IVa__IIa_IVa) { IIb_IVa } else { IIIa_IVb }
            
        }
    }

    ## Compute parameters for hat and squeeze for initial intervals.
    for (i in 1:(n.ivs+1)) {
        if (is.finite(params["next",i])) {
            params[,i] <- hat.iv(left=params[,i], right=params[,params["next",i]])
        }
    }
    
    ## DEBUG:
    ## print(params[,1:(n.ivs+1)])   ## DEBUG

    ## Compute total areas for initial hat and squeeze.
    A.ht.tot <- sum(params["A.ht", 1:n.ivs], na.rm=TRUE)
    A.sq.tot <- sum(params["A.sq", 1:n.ivs], na.rm=TRUE)

    ## We have to split intervals where the area
    ## between hat and squueze is too large.
    while (! isTRUE(A.ht.tot / A.sq.tot <= rho)) {
        
        ## Compute average area.
        threshold <- 0.99 * (A.ht.tot - A.sq.tot) / n.ivs

        ## Find all intervals where area between hat and squeeze is too large.
        toolarge <- which(params["A.ht", 1:(n.ivs+1)] - params["A.sq", 1:(n.ivs+1)] >= threshold)
        if (length(toolarge) <= 0) break;
        if (n.ivs + length(toolarge) >= max.intervals) {
            warning ("maximum number of intervals exceeded")
            break;
        }
    
        ## Split all intervals.
        for (i in toolarge) {
            ## Get boundary points of interval.
            left <- params[,i]
            right <- params[,params["next",i]]

            ## Splitting point and types of subinervals
            split.iv <- type.no2ndD.split(left=left, right=right, lpdf, dlpdf)
            p <- split.iv[1]
            Tfp <- split.iv[2]
            dTfp <- split.iv[3]
            typel <- split.iv[4]
            typer <- split.iv[5]

            ## Increment number of intervals
            n.ivs <- n.ivs+1
            
            ## Create new interval on r.h.s.
            params["x", n.ivs+1] <- p
            params["c", n.ivs+1] <- left["c"]
            params[c("Tfx","dTfx","d2Tfx"), n.ivs+1] <- c(Tfp, dTfp, NA_real_)
            ## or equivalently:
            ##   params[c("Tfx","dTfx","d2Tfx"), n.ivs+1] <- Tfdd(lpdf, dlpdf, NULL, left["c"], p)
            params["type",n.ivs+1] <- typer
            params[      ,n.ivs+1] <- hat.iv(left=params[,n.ivs+1], right=right)
            params["next",n.ivs+1] <- left["next"]

            ## Update interval on l.h.s.
            params["type",i] <- typel
            params[      ,i] <- hat.iv(left=params[,i], right=params[,n.ivs+1])
            params["next",i] <- n.ivs+1
        }
       
        ## Update total areas.
        A.ht.tot <- sum(params["A.ht",], na.rm=TRUE)
        A.sq.tot <- sum(params["A.sq",], na.rm=TRUE)
    }

    ## DEBUG: 
    ## print(params[,1:(n.ivs+1)])

    ## Check result.
    if (! is.finite(A.ht.tot)) {
        stop ("Cannot create hat function. A_hat is not finite: A_hat=", A.ht.tot)
    } 
    
    if (! isTRUE(A.ht.tot > 0.)) {
        stop ("Cannot create hat function. A_hat is not > 0: A_hat=", A.ht.tot)
    } 
    
    if (! is.finite(A.sq.tot)) {
        stop ("Cannot create squeeze. A_squeeze is not finite: A_squeeze=", A.sq.tot)
    } 
    
    if (! isTRUE(A.sq.tot >= 0.)) {
        stop ("Cannot create squeeze. A_squeeze is not >= 0: A_squeeze=", A.sq.tot)
    } 
    
    if (isTRUE(A.ht.tot < A.sq.tot))
        stop ("Invalid input: A_hat < A_squeeze!")
    
    if (!isTRUE(A.ht.tot / A.sq.tot <= rho))
        warning ("ratio A_hat / A_squeeze = ", A.ht.tot / A.sq.tot,
                 " is larger than rho = ",rho)
    
    ## Truncate working array.
    params <- params[,1:(n.ivs+1)]
    
    ## Sort array (so that the linked list becomes obsolete).
    params <- params[,order(params["x",])]
    
    ## Compute guide table.
    ## For R version use
    ##   guidetable <- make.guidetable(params)
    ##   Acum <- guidetable$Acum
    ##   gt <- guidetable$gt
    ## C version:
    Acum <- numeric(n.ivs)
    gt <- integer(n.ivs)
    A.ht.tot <- .Call("Tinflex_RC_make_guide_table", params, Acum, gt)
    
    ## Create S3 class that contains generator.
    generator <- list(
        ivs=params,          ## data for hat and squeeze
        lpdf=lpdf,           ## log-density
        A.ht.tot=A.ht.tot,   ## total area below hat
        A.sq.tot=A.sq.tot,   ## total area below hat
        env=parent.frame(2), ## environment for evaluating log-density
        iniv=iniv,           ## initial intervals (for print.Tinflex)
        Acum=Acum,           ## cumulated areas
        gt=gt,               ## guide table
        haved2Tf=FALSE       ## whether 2nd derivative is used (for print.Tinflex)
    )
    class(generator) <- "Tinflex"
    
    ## Return generator object. 
    return (generator)
} ## -- end of Tinflex.setup.no2ndD() -- ##

## --------------------------------------------------------------------------

