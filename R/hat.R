#############################################################################
##
##  Parameters for hat and squeeze in an interval.
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
##  If 'y' is either the left boundary or the right boundary.
##  If possible then 'y' is chosen as the boundary point where the transformed
##  obtains higher valuation.

## --------------------------------------------------------------------------

hat.iv <- function(left, right) {
    ## ----------------------------------------------------------------------
    ## Compute hat and squeeze for a paricular interval.
    ## ----------------------------------------------------------------------
    ##   left   ... data for boundary point to the left and for entire interval
    ##   right  ... data for boundary point to the right
    ## ----------------------------------------------------------------------
    ## Return: vector with parameter of same kind as 'left'.
    ## ----------------------------------------------------------------------

    ## Check for interval of length 0.
    if (isTRUE( left["x"] == right["x"] )) {
        left[c("ht.a","ht.b","ht.y")] <- c(left["Tfx"], 0, left["x"])
        left[c("sq.a","sq.b","sq.y")] <- c(left["Tfx"], 0, left["x"])
        left["A.ht"] <- 0
        left["A.sq"] <- 0
        left["type"] <- 0
        return (left)
    }
    
    ## Get type of distribution.
    type <- left["type"]

    ## Unknown type
    if (type == 0) {
        left["A.ht"] <- Inf
        left["A.sq"] <- 0
        return (left)
    }
    
    ## Compute tangent at boundary points.
    tl <- c( left["Tfx"], left["dTfx"], left["x"])
    tr <- c( right["Tfx"], right["dTfx"], right["x"])
    
    ## Compute secant.
    R <- (right["Tfx"] - left["Tfx"]) / (right["x"] - left["x"])
    if (isTRUE(left["Tfx"] >= right["Tfx"])) {
        sc <- c( left["Tfx"], R, left["x"])
    }
    else {
        sc <- c( right["Tfx"], R, right["x"])
    }
    
    ## Case: unbounded domains.
    if (is.infinite(left["x"]) && type == IVa) {
        left[c("ht.a","ht.b","ht.y")] <- tr
        left[c("sq.a","sq.b","sq.y")] <- NA
    }
    else if (is.infinite(right["x"]) && type == IVa) {
        left[c("ht.a","ht.b","ht.y")] <- tl
        left[c("sq.a","sq.b","sq.y")] <- NA
    }
    
    ## Case: bounded domains.
    else if (type == Ia) {
        left[c("ht.a","ht.b","ht.y")] <- tl
        left[c("sq.a","sq.b","sq.y")] <- tr
    }
    else if (type == Ib) {
        left[c("ht.a","ht.b","ht.y")] <- tr
        left[c("sq.a","sq.b","sq.y")] <- tl
    }
    else if (type == IIa) {
        left[c("ht.a","ht.b","ht.y")] <- tl
        left[c("sq.a","sq.b","sq.y")] <- sc
    }
    else if (type == IIb) {
        left[c("ht.a","ht.b","ht.y")] <- tr
        left[c("sq.a","sq.b","sq.y")] <- sc
    }
    else if (type == IIIa) {
        left[c("ht.a","ht.b","ht.y")] <- sc
        left[c("sq.a","sq.b","sq.y")] <- tr
    }
    else if (type == IIIb) {
        left[c("ht.a","ht.b","ht.y")] <- sc
        left[c("sq.a","sq.b","sq.y")] <- tl
    }
    else if (type == IVa) {
        if (isTRUE(left["Tfx"] > right["Tfx"]))
            left[c("ht.a","ht.b","ht.y")] <- tl
        else
            left[c("ht.a","ht.b","ht.y")] <- tr
        left[c("sq.a","sq.b","sq.y")] <- sc
    }
    else if (type == IVa) {
        if (isTRUE(left["Tfx"] > right["Tfx"]))
            left[c("ht.a","ht.b","ht.y")] <- tl
        else
            left[c("ht.a","ht.b","ht.y")] <- tr
        left[c("sq.a","sq.b","sq.y")] <- sc
    }
    else if (type == IVb) {
        left[c("ht.a","ht.b","ht.y")] <- sc
        if (isTRUE(left["Tfx"] > right["Tfx"]))
            left[c("sq.a","sq.b","sq.y")] <- tr
        else
            left[c("sq.a","sq.b","sq.y")] <- tl
    }

    ## Case: combined types
    else if (type == IIa_IVa) {
        left[c("ht.a","ht.b","ht.y")] <- tl
        left[c("sq.a","sq.b","sq.y")] <- sc
    }
    else if (type == IIb_IVa) {
        left[c("ht.a","ht.b","ht.y")] <- tr
        left[c("sq.a","sq.b","sq.y")] <- sc
    }
    else if (type == IIIa_IVb) {
        left[c("ht.a","ht.b","ht.y")] <- sc
        left[c("sq.a","sq.b","sq.y")] <- tr
    }
    else if (type == IIIb_IVb) {
        left[c("ht.a","ht.b","ht.y")] <- sc
        left[c("sq.a","sq.b","sq.y")] <- tl
    }

    else {  ## Something went wrong
        left["A.ht"] <- Inf
        left["A.sq"] <- 0
        return (left)
    }

    ## Compute area below hat.
    left["A.ht"] <- area(left["c"], left["ht.a"], left["ht.b"], left["ht.y"], left["x"], right["x"])
    
    ## Compute area below squeeze.
    A.sq <- area(left["c"], left["sq.a"], left["sq.b"], left["sq.y"], left["x"], right["x"])
    left["A.sq"] <- if (is.finite(A.sq)) { A.sq } else { 0 } 

    ## Return vector with parameters.
    return (left)

}  ## -- end of hat.iv() -- ##

## --------------------------------------------------------------------------
