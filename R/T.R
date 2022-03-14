#############################################################################
##
##  Transformed Density.
##
#############################################################################

Tfdd <- function(lpdf, dlpdf, d2lpdf, cT, x) {
    ## ----------------------------------------------------------------------
    ## Compute transformed density and its derivatives by means of
    ## log-density and its derivatives.
    ## If d2lpdf == NULL then the 2nd derivative is set to NA
    ## ----------------------------------------------------------------------
    ##   lpdf   ... log-density 
    ##   dlpdf  ... derivative of log-density
    ##   d2lpdf ... 2nd derivative of log-density
    ##   cT     ... parameter for transformation
    ##   x      ... argument
    ## ----------------------------------------------------------------------
    ## Return: vector consisting of
    ##   Tfx    ... transformed density at x
    ##   dTfx   ... derivative of transformed density at x
    ##   d2Tfx  ... second derivative of transformed density at x
    ## ----------------------------------------------------------------------
    
    ## Evaluate log density and its derivatives.
    lfx <- lpdf(x)
    dlfx <- dlpdf(x)
    d2lfx <- if (is.null(d2lpdf)) { NA_real_ } else { d2lpdf(x) }
    
    if (cT == 0) {
        ## Case: T(x) = log(x)
        return (c(lfx,dlfx,d2lfx))
    }
    else {
        ## Case: T(x) = sign(c) * x^c
        Tfx   <- sign(cT) * exp(cT * lfx)
        dTfx  <- cT * Tfx * dlfx
        d2Tfx <- if (is.null(d2lpdf)) { NA_real_ } else { cT * Tfx * (cT * dlfx^2 + d2lfx) }
        return (c(Tfx,dTfx,d2Tfx))
    }

}  ## -- end of Tfdd() -- ## 

## --------------------------------------------------------------------------

Tf <- function(lpdf, cT, x) {
    ## ----------------------------------------------------------------------
    ## Compute transformed density by means of its log-density.
    ## ----------------------------------------------------------------------
    ##   lpdf   ... log-density 
    ##   cT     ... parameter for transformation
    ##   x      ... argument
    ## ----------------------------------------------------------------------
    ## Return: transformed density at x.
    ## ----------------------------------------------------------------------
    
    if (cT == 0) { lpdf(x) } else { sign(cT) * exp(cT * lpdf(x)) }

}  ## -- end of Tf() -- ##

## --------------------------------------------------------------------------

Tinv <- function(cT, x) {
    ## ----------------------------------------------------------------------
    ## Compute inverse transformation.
    ## ----------------------------------------------------------------------
    ##   cT     ... parameter for transformation
    ##   x      ... argument
    ## ----------------------------------------------------------------------
    ## Return: inverse transformation at x.
    ## ----------------------------------------------------------------------
    
    if (cT == 0) {
        ## Case: T(x) = log(x)
        return (exp(x))
    }
    else if (cT == -0.5) {
        ## Case: T(x) = -1/sqrt(x)
        return (1/(x*x))
    }
    else if (cT == 1) {
        ## Case: T(x) = x
        return (x)
    }
    else {
        ## Case: T(x) = sign(c) * x^c
        return ((sign(cT)*x)^(1/cT))
    }
    
}  ## -- end of Tinv() -- ##

## --------------------------------------------------------------------------

FT <- function(cT, x) {
    ## ----------------------------------------------------------------------
    ## Compute antiderivative of inverse transformation.
    ## ----------------------------------------------------------------------
    ##   cT     ... parameter for transformation
    ##   x      ... argument
    ## ----------------------------------------------------------------------
    ## Return: antiderivative of inverse transformation at x.
    ## ----------------------------------------------------------------------
    
    if (cT == 0) {
        ## Case: T(x) = log(x)
        return (exp(x))
    }
    else if (cT == -0.5) {
        ## Case: T(x) = -1/sqrt(x)
        return (-1/x)
    }
    else if (cT == -1) {
        ## Case: T(x) = -1/x
        return (-log(-x))
    }
    else {
        ## Case: T(x) = sign(c) * x^c
        return (sign(cT) * cT/(cT+1) * (sign(cT) * x)^((cT+1)/cT))
    }
    
}  ## -- end of FT() -- ##

## --------------------------------------------------------------------------

FTinv <- function(cT, x) {
    ## ----------------------------------------------------------------------
    ## Compute inverse of antiderivative of inverse transformation.
    ## ----------------------------------------------------------------------
    ##   cT     ... parameter for transformation
    ##   x      ... argument
    ## ----------------------------------------------------------------------
    ## Return: inverse of antiderivative of inverse transformation at x.
    ## ----------------------------------------------------------------------
    
    if (cT == 0) {
        ## Case: T(x) = log(x)
        return (log(x))
    }
    else if (cT == -0.5) {
        ## Case: T(x) = -1/sqrt(x)
        return (-1/x)
    }
    else if (cT == -1) {
        ## Case: T(x) = -1/x
        return (-exp(-x))
    }
    else {
        ## Case: T(x) = sign(c) * x^c
        return (sign(cT) * (sign(cT) * (cT+1)/cT * x)^(cT/(cT+1)))
    }

}  ## -- end of FTinv() --##

## --------------------------------------------------------------------------
