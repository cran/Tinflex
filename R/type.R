#############################################################################
##
##  Types of intervals.
##
#############################################################################

## --------------------------------------------------------------------------
## Codes for types of interval (see paper):

Ia   <- -1
Ib   <-  1
IIb  <-  2
IIa  <- -2
IIIa <- -3
IIIb <-  3
IVa  <- -4    ## concave
IVb  <-  4    ## convex

## Codes for combined types (see new paper):

IIa_IVa  <-  -24   ## -2 -4
IIb_IVa  <-   24   ##  2 -4
IIIa_IVb <-  -34   ## -3  4 
IIIb_IVb <-   34   ##  3  4

## split types (see new paper):

IIb_IVa__IIa_IVa   <-  222    ## For case (3.3.3)
IIIa_IVb__IIIb_IVb <-  333    ## For case (4.3.3)


## --------------------------------------------------------------------------

type.rules <- function(dTfl,dTfr, d2Tfl,d2Tfr, R) {
    ## ----------------------------------------------------------------------
    ## Estimate type of interval.
    ## ----------------------------------------------------------------------
    ##   dTfl   ... derivative on left boundary
    ##   dTfr   ... derivative on right boundary
    ##   d2Tfl  ... 2nd derivative on left boundary
    ##   d2Tfr  ... 2nd derivative on right boundary
    ##   R      ... slope of secant
    ## ----------------------------------------------------------------------
    ## Return: type of interval;
    ##         0 if type cannot be estimated.
    ## ----------------------------------------------------------------------

    ## Regular types
    
    if (isTRUE(dTfl > R && dTfr > R)) {
        return (Ia)
    }

    if (isTRUE(dTfl < R && dTfr < R)) {
        return (Ib)
    }

    if (isTRUE(d2Tfl < 0 && d2Tfr < 0)) {
        return (IVa)
    }
    
    if (isTRUE(d2Tfl > 0 && d2Tfr > 0)) {
        return (IVb)
    }

    if (isTRUE(d2Tfl <= 0 && 0 <= d2Tfr)) {
        if (isTRUE(dTfl >= R && R >= dTfr)) {
            return (IIa)
        }
        if (isTRUE(dTfl <= R && R <= dTfr)) {
            return (IIIa)
        }
    }

    if (isTRUE(d2Tfl >= 0 && 0 >= d2Tfr)) {
        if (isTRUE(dTfl >= R && R >= dTfr)) {
            return (IIb)
        }
        if (isTRUE(dTfl <= R && R <= dTfr)) {
            return (IIIb)
        }
    }

    ## Combinded types

    if (isTRUE(dTfl >= R && R >= dTfr)) {
        if (isTRUE(d2Tfl <= 0)) {
            return (IIa_IVa)
        }
        if (isTRUE(d2Tfr <= 0)) {
            return (IIb_IVa)
        }
        ## else: should not happen
    }
    
    if (isTRUE(dTfl <= R && R <= dTfr)) {
        if (isTRUE(d2Tfr >= 0)) {
            return (IIIa_IVb)
        }
        if (isTRUE(d2Tfl >= 0)) {
            return(IIIb_IVb)
        }
        ## else: should not happen
    }

    ## Cannot estimate type of interval.
    ## Do we need a warning? Probably not.
    ## --> warning("cannot detect type of interval")

    return (0)

} ## -- end of type.rules() -- ##

## --------------------------------------------------------------------------

type.with2ndD <- function(left, right) {
    ## ----------------------------------------------------------------------
    ## Estimate type of interval.
    ## ----------------------------------------------------------------------
    ##   left   ... data for boundary point to the left
    ##   right  ... data for boundary point to the right
    ## ----------------------------------------------------------------------
    ## Return: type of interval;
    ##         0 if type cannot be estimated.
    ## ----------------------------------------------------------------------
    
    ## Case: Unbounded domain.
    ## Assumption: transformed density is concave
    ## (We check at the finite boundary)

    if ( is.infinite(left["x"])) {
        if (isTRUE(right["d2Tfx"] < 0 && right["dTfx"] >= 0)) {
            return (IVa)
        }
        else {
            return (0)
        }
    }
    
    if ( is.infinite(right["x"])) {
        if (isTRUE(left["d2Tfx"] < 0 && left["dTfx"] <= 0)) {
            return (IVa)
        }
        else {
            return (0)
        }
    }

    ## Parameter 'c' for interval.
    cT <- left["c"]
        
    ## Case: Interval where the density vanishes at boundary
    ##       (and thus the log-density is -Inf).
    ## Subcase: c <= 0

    if (isTRUE(cT<=0 && left["Tfx"] == -Inf)) {
        if (isTRUE(right["d2Tfx"] < 0 && right["dTfx"] >= 0)) {
            return (IVa)
        }
        else {
            return (0)
        }
    }

    if (isTRUE(cT<=0 && right["Tfx"] == -Inf)) {
        if (isTRUE(left["d2Tfx"] < 0 && left["dTfx"] <= 0)) {
            return (IVa)
        }
        else {
            return (0)
        }
    }

    ## Subcase: c > 0
    ## We have no information about the derivative at the boundary point,
    ## as we use the log-density as argument and not the density itself.
    ## Assumption: Tf does not have an inflection point in this interval and
    ## is thus either concave or convex in this interval.
    if (isTRUE(cT>0 && left["Tfx"] == 0)) {
        if (isTRUE(right["d2Tfx"] < 0 && right["dTfx"] >= 0)) {
            return (IVa)
        }
        if (isTRUE(right["d2Tfx"] > 0 && right["dTfx"] >= 0)) {
            return (IVb)
        }
        else {
            return (0)
        }
    }
    
    if (isTRUE(cT>0 && right["Tfx"] == 0)) {
        if (isTRUE(left["d2Tfx"] < 0 && left["dTfx"] <= 0)) {
            return (IVa)
        }
        if (isTRUE(left["d2Tfx"] > 0 && left["dTfx"] <= 0)) {
            return (IVb)
        }
        else {
            return (0)
        }
    }

    ## Case: Domains where the density has a pole at boundary
    ##       (and thus the transformed density equals 0 for c<0).
    ## Assumption: the transformed density is convex.
    ## (We make a simple check)
    if (isTRUE(cT<0)) {
        if ((left["Tfx"] == 0  && right["d2Tfx"] > 0) ||
            (right["Tfx"] == 0 && left["d2Tfx"] > 0)  ) {
            return (IVb)
        }
    }

    ## General case
    
    ## Compute slope of secant.
    R <- (right["Tfx"] - left["Tfx"]) / (right["x"] - left["x"])

    type <- type.rules(left["dTfx"], right["dTfx"],
               left["d2Tfx"], right["d2Tfx"], R)

    if (! (abs(type) < 10)) {
        ## we get a combined type. This should not happen.
        ## So we simply return "not determined".
        return(0)
    }

    return(type)
    
}  ## -- end of type.with2ndD() -- ## 

## --------------------------------------------------------------------------

type.no2ndD.init <- function(left, right, lpdf, dlpdf) {
    ## ----------------------------------------------------------------------
    ## Estimate type of interval from scratch.
    ## ----------------------------------------------------------------------
    ##   left   ... data for boundary point to the left
    ##   right  ... data for boundary point to the right
    ##   lpdf   ... log-density 
    ##   dlpdf  ... derivative of log-density
    ## ----------------------------------------------------------------------
    ## Return: type of interval;
    ##         0 if type cannot be estimated.
    ## ----------------------------------------------------------------------

    ## Case: interval of length 0.
    if (isTRUE(left["x"] == right["x"] )) {
        return (0)
    }

    ## Parameter 'c' for interval.
    cT <- left["c"]

    ## Control point
    p <- arc.mean(left["x"],right["x"])
    tmp <- Tfdd(lpdf, dlpdf, NULL, cT, p)
    Tfp <- tmp[1]
    dTfp <- tmp[2]
    
    ## Case: Unbounded domain.
    ## Assumption: transformed density is concave
    ## (We make a simple check for this condition.)

    if ( is.infinite(left["x"])) {
        if (isTRUE(dTfp >= right["dTfx"] && right["dTfx"] > 0)) {
            return (IVa)
        } else {
            return (IIa_IVa)
        }
    }
    
    if ( is.infinite(right["x"])) {
        if (isTRUE(left["dTfx"] >= dTfp && left["dTfx"] < 0)) {
            return (IVa)
        } else {
            return (IIb_IVa)
        }
    }

    ## Evaluate tangents at boundary points at p
    tlp <- as.numeric(left["Tfx"] + left["dTfx"] * (p - left["x"]))
    trp <- as.numeric(right["Tfx"] + right["dTfx"] * (p - right["x"]))

    ## Compute slope of secant.
    R <- (right["Tfx"] - left["Tfx"]) / (right["x"] - left["x"])


    ## Case: Interval where the density vanishes at boundary
    ##       (and thus the log-density is -Inf).
    ## Subcase: c <= 0

    if (isTRUE(cT <= 0 && left["Tfx"] == -Inf)) {
        ## print("case (Inf-IIa_IVa)")
        return (IIa_IVa)
    }
    if (isTRUE(cT <= 0 && right["Tfx"] == -Inf)) {
        ## print("case (Inf-IIb_IVa)")
        return (IIb_IVa)
    }

    ## Subcase: c > 0
    ## We have no information about the derivative at the boundary point,
    ## as we use the log-density as argument and not the density itself.
    ## Assumption: Tf does not have an inflection point in this interval and
    ## is thus either concave or convex in this interval.
    if (isTRUE(cT > 0 && left["Tfx"] == 0 && right["Tfx"] > 0)) {
        if (isTRUE(right["dTfx"] <= R)) {
            ## concave
            return (IVa);
        }
        else if (isTRUE(right["dTfx"] >= R)) {
            ## convex
            return (IVb);
        }
        else {
            return (0)
        }
    }

    if (isTRUE(cT > 0 && left["Tfx"] > 0 && right["Tfx"] == 0)) {
        if (isTRUE(left["dTfx"] >= R)) {
            ## concave
            return (IVa)
        }
        else if (isTRUE(left["dTfx"] <= R)) {
            ## convex
            return (IVb)
        }
        else {
            return (0)
        }
    }

    ## Case: Domains where the density has a pole at boundary
    ##       (and thus the transformed density equals 0 for c<0).
    ## Assumption: the transformed density is convex.
    ## (We make a simple check)
    if (isTRUE(cT<0)) {
        if ((isTRUE(left["Tfx"] == 0 && R <= right["dTfx"] && right["dTfx"] < 0)) ||
            (isTRUE(right["Tfx"] == 0 && R >= left["dTfx"] && left["dTfx"] > 0)) ) {
            return (IVb)
        }
    }
    
    ## Check for all other possible cases.
    if (isTRUE(left["dTfx"] >= R && right["dTfx"] >= R)) {
        ## print("case (1)")
        return (Ia)
    }
    
    if (isTRUE(left["dTfx"] <= R && right["dTfx"] <= R)) {
        ## print("case (2)")
        return (Ib)
    }

    if (isTRUE(left["dTfx"] >= R && right["dTfx"] <= R)) {
        ## print("case (3)")
        if (isTRUE(dTfp <= right["dTfx"])) {
            ## print("case (3.1)")
            return (IIa)
        }
        if (isTRUE(dTfp >= left["dTfx"])) {
            ## print("case (3.2)")
            return (IIb)
        }
        ## else: left["dTfx"] >= dTfp >= right["dTfx"]
        ## case (3.3)
        if (isTRUE(Tfp > tlp)) {
            ## print("case (3.3.1)")
            return (IIb)
        }
        if (isTRUE(Tfp > trp)) {
            ## print("case (3.3.2)")
            return (IIa)
        }
        else {
            ## print("case (3.3.3)")
            if (! isTRUE(Tfp <= tlp && Tfp <= trp)) {
                ## check for possible numerical errors
                return (0)
            }
            return (IIb_IVa__IIa_IVa)
        }
    }

    if (isTRUE(left["dTfx"] <= R && right["dTfx"] >= R)) {
        ## print("case (4)")
        if (isTRUE(dTfp <= left["dTfx"])) {
            ## print("case (4.1)")
            return (IIIa)
        }
        if (isTRUE(dTfp >= right["dTfx"])) {
            ## print("case (4.2)")
            return (IIIb)
        }
        ## else: left["dTfx"] <= dTfp <= right["dTfx"]
        ## case (4.3)
        if (isTRUE(Tfp < tlp)) {
            ## print("case (4.3.1)")
            return (IIIa)
        }
        if (isTRUE(Tfp < trp)) {
            ## print("case (4.3.2)")
            return (IIIb)
        }
        else {
            ## print("case (4.3.3)")
            if (! isTRUE(Tfp >= tlp && Tfp >= trp)) {
                ## check for possible numerical errors
                return (0)
            }
            return (IIIa_IVb__IIIb_IVb)
        }
    }

    ## Cannot estimate type of interval.
    ## Do we need a warning? Probably not.
    ## > warning("cannot detect type of interval")
    
    return (0)

} ## -- end of type.no2ndD.init() -- ##

## ..........................................................................

type.no2ndD.split <- function(left, right, lpdf, dlpdf) {
    ## ----------------------------------------------------------------------
    ## Estimate types of splitted subintervals.
    ## ----------------------------------------------------------------------
    ##   left   ... data for boundary point to the left
    ##   right  ... data for boundary point to the right
    ##   lpdf   ... log-density 
    ##   dlpdf  ... derivative of log-density
    ## ----------------------------------------------------------------------
    ## Return:
    ##   array: splitting point, Tf(p), dTf(p),
    ##          type left interval, type right interval
    ## ----------------------------------------------------------------------

    ## Type of interval to be split.
    type <- left["type"]

    ## Splitting point
    p <- arc.mean(left["x"],right["x"])
    tmp <- Tfdd(lpdf, dlpdf, NULL, left["c"], p)
    Tfp <- tmp[1]
    dTfp <- tmp[2]

    ## Case: Unknown type in current interval
    if (type == 0) {
        ## try initial interval
        type <- type.no2ndD.init(left, right, lpdf, dlpdf)

        if (type == IIb_IVa__IIa_IVa) {           ## Case (3.3.3)
            return(c(p, Tfp, dTfp, IIb_IVa, IIa_IVa))
        }
        else if (type == IIIa_IVb__IIIb_IVb ) {   ## Case (4.3.3)
            return(c(p, Tfp, dTfp, IIIa_IVb, IIIb_IVb))
        }
        ## else: continue with new type
    }

    ## Concave | convex
    if (type == IVa) {
        return(c(p, Tfp, dTfp, IVa, IVa))
    }

    if (type == IVb) {
        return(c(p, Tfp, dTfp, IVb, IVb))
    }

    ## Auxiliary splitting point
    if (is.finite(left["x"]) && is.finite(right["x"])) {
        pD <- p + (right["x"] - left["x"]) * 1e-3
    }
    else if (is.finite(left["x"])) {
        pD <- p + (p - left["x"]) * 1e-3
    }
    else {  ## is.finite(right["x"])
        pD <- p + (right["x"] - p) * 1e-3
    }
    tmp <- Tfdd(lpdf, dlpdf, NULL, left["c"], pD)
    TfpD <- tmp[1]
    dTfpD <- tmp[2]

    if (! all(is.finite(c(pD, TfpD, dTfpD)))) {
        ## something is wrong here
        return (c(p,Tfp,dTfp,0,0))
    }

    ## Get splitting points and sign of 2nd derivative

    psplit <- NA
    d2Tfl <- NA
    d2Tfr <- NA

    if (type == Ia ||
        type == IIa || type == IIIa) {
        d2Tfl <- -1
        d2Tfr <- +1
        if (dTfp < dTfpD) {
            psplit <- "pD"
            d2Tfp <- +1
        }
        else {
            psplit <- "p"
            d2Tfp <- -1
        }
    }

    else if (type == Ib || type == IIb || type == IIIb) {
        d2Tfl <- +1
        d2Tfr <- -1
        if (dTfp < dTfpD) {
            psplit <- "p"
            d2Tfp <- +1
        }
        else {
            psplit <- "pD"
            d2Tfp <- -1
        }
    }

    ## Combined types
    
    else if (type == IIa_IVa) {
        d2Tfl <- -1
        d2Tfr <- NA
        if (dTfp < dTfpD) {
            psplit <- "pD"
            d2Tfp <- +1
            d2Tfr <- +1
        }
        else {
            psplit <- "p"
            d2Tfp <- -1
        }
    }

    else if (type == IIb_IVa) {
        d2Tfl <- NA
        d2Tfr <- -1
        if (dTfp < dTfpD) {
            psplit <- "p"
            d2Tfp <- +1
            d2Tfl <- +1
        }
        else {
            psplit <- "pD"
            d2Tfp <- -1
        }
    }

    else if (type == IIIa_IVb) {
        d2Tfl <- NA
        d2Tfr <- +1
        if (dTfp < dTfpD) {
            psplit <- "pD"
            d2Tfp <- +1
        }
        else {
            psplit <- "p"
            d2Tfl <- -1
            d2Tfp <- -1
        }
    }

    else if (type == IIIb_IVb) {
        d2Tfl <- +1
        d2Tfr <- NA
        if (dTfp < dTfpD) {
            psplit <- "p"
            d2Tfp <- +1
        }
        else {
            psplit <- "pD"
            d2Tfp <- -1
            d2Tfr <- -1
        }
    }

    else {
        ## this should not happen
        return (c(p,Tfp,dTfp,0,0))
    }

    
    ## If you have to use splitting point 'pD'
    ## we copy it into 'p' first
    if (psplit != "p") {
        p <- pD
        Tfp <- TfpD
        dTfp <- dTfpD
    }
    
    ## Compute slopes of secants in the subintervals.
    Rl <- (Tfp - left["Tfx"]) / (p - left["x"])
    Rr <- (right["Tfx"] - Tfp) / (right["x"] - p)

    ## Get types for each subinterval
    typel <- type.rules(left["dTfx"], dTfp, d2Tfl, d2Tfp, Rl) 
    typer <- type.rules(dTfp, right["dTfx"], d2Tfp, d2Tfr, Rr) 

    ## return result
    c(p, Tfp, dTfp, typel, typer)

} ## -- end of type.no2ndD.split() -- ##

## --------------------------------------------------------------------------
