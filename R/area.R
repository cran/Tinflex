#############################################################################
##
##  Area below hat and squeeze
##
#############################################################################

##  Remark:
##  Tangents and secants have the form
##
##    t(x) = a + b*(x-y)
##
##  where 
##    a ... intercept for hat in transformed scale
##    b ... slope for hat in transformed scale
##    y ... anchor point of linear function
##

## --------------------------------------------------------------------------

area <- function(cT, a,b,y, from,to) {
    ## ----------------------------------------------------------------------
    ## Compute area underneath hat or squeeze in particular interval.
    ## ----------------------------------------------------------------------
    ##   cT       ... parameter for transformation
    ##   a, b, y  ... intercept, slope and anchor point of transformed hat or squeeze
    ##   from, to ... range of integration
    ## ----------------------------------------------------------------------
    ## Return: area;
    ##         0 in case of an error.
    ## ----------------------------------------------------------------------

    ## Compute area ...
    area <- do.area(cT,a,b,y,from,to)
    
    ## ... and check result.
    
    if (! is.finite(area)) {
        ## We return Inf in all cases where 'area' is also not finite (e.g. NaN or NA)
        area <- Inf
    }
    
    if (area < 0) {
        ## Area is strictly negative.
        ## Two possible reasons:
        ## 1. The correct value of 'area' is extremely small or 0.
        ##    The negative number might result from cancelation.
        ##    Then 'area' can be set to 0.
        ## 2. Due to sever round-off error, the result is numerical garbage.
        ##
        ## We observed (2) in our experiments. So we discard the result and
        ## return 'Inf' (which means that the interval must be split).
        
        area <- Inf
    }
    
    return (area)

}  ## -- end of area() -- ##

## ..........................................................................

do.area <- function(cT, a,b,y, from,to) {
    ## ----------------------------------------------------------------------
    ## Perform computation of area underneath hat or squeeze.
    ## ----------------------------------------------------------------------
    ##   cT       ... parameter for transformation
    ##   a, b, y  ... intercept, slope and anchor point of transformed hat or squeeze
    ##   from, to ... range of integration
    ## ----------------------------------------------------------------------
    ## Return: area.
    ## ----------------------------------------------------------------------

    ## check for a new "empty" interval without data
    if (is.na(a) && !is.nan(a)) {
        ## We have a new interval without data.
        return (Inf)
    }

    ## Test where the tangent is constructed:
    ##   s = +1 if tangent is constructed on lower boundary of interval;
    ##       -1 if tangent is constructed on upper boundary of interval.
    s <- if (isTRUE((to-y)>(y-from))) 1 else -1
    
    ## Generally we have
    ##   area <- (FT(cT, a+b*(to-y)) - FT(cT, a+b*(from-y))) / b
    ## For numerical reasons we have to distinguish
    ## between different values of 'cT'.
    
    if (cT == 0) {
        ## Case: T(x)=log(x)
        
        z <- s * b*(to-from)
        if (isTRUE(abs(z) > 1.e-6)) {
            area <- (exp(a+b*(to-y)) - exp(a+b*(from-y))) / b
        } else {
            ## We need approximation by Taylor polynomial to avoid
            ## severe round-off errors.
            area <- exp(a) * (to-from) * (1 + z/2 + z*z/6)
        }
        return (area)
    }
    
    ## else: c!=0
    
    ## The tangent to the transformed density must not vanish.
    ## Otherwise, we cannot construct the hat function.
    tleft  <- a+b*(from-y)
    tright <- a+b*(to-y)
    
    if (isTRUE(cT<0) && !isTRUE(tleft <= 0 && tright <= 0) ) {
        ## we have to split the interval
        return (Inf)
    }
    else if (isTRUE(cT>0) && !isTRUE(tleft >= 0 && tright >= 0) ) {
        return (Inf)
    }

    ## Transform b.
    z <- s * b/a * (to-from)
    
    if (cT == -0.5) {
        ## Case: T(x) = -1/sqrt(x)
        
        if (isTRUE(abs(z) > 0.5)) {
            area <- (-1/(a+b*(to-y)) + 1/(a+b*(from-y))) / b
        } else {
            area <- 1/(a*a) * (to-from) / (1 + z)
        }
        return (area)
    }
    
    if (cT == -1) {
        ## Case: T(x) = -1/x
        
        if (isTRUE(abs(z) > 1.e-6)) {
            area <- (-log(-a-b*(to-y)) + log(-a-b*(from-y))) / b
        } else {
            ## Approximation by Taylor polynomial.
            area <- -1/a * (to-from) * (1 - z/2 + z*z/3)
        }
        return (area)
    }
    
    if (cT == 1) {
        ## Case: T(x) = x
        
        area <- 0.5 * a * (to-from) * (z+2)
        return (area)
    }
    
    ## Case: T(x) = sgn(c) * x^c
    ## For all other cases we only use a rough approximation in
    ## case of numerical errors.
    
    if (isTRUE(abs(b)>1e-10)) {
        area <- (FT(cT, a+b*(to-y)) - FT(cT, a+b*(from-y))) / b
    } else {
        area <- Tinv(cT, a) * (to-from)
    }
    return (area)
    
}  ## -- end of do.area() -- ##

## --------------------------------------------------------------------------
