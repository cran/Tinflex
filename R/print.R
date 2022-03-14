#############################################################################
##
##  Print data about generator.
##
#############################################################################

print.Tinflex <- function(x, debug=FALSE, ...) {
    ## ----------------------------------------------------------------------
    ## S3 method for printing class 'Tinflex'.
    ## ----------------------------------------------------------------------
    ##   x     ... S3 object of class 'Tinflex'
    ##   debug ... if TRUE print intervals
    ## ----------------------------------------------------------------------
    
    Tinflex.print(x, debug, classname="Tinflex", ...)
}

## --------------------------------------------------------------------------

print.TinflexC <- function(x, debug=FALSE, ...) {
    ## ----------------------------------------------------------------------
    ## S3 method for printing class 'TinflexC'.
    ## ----------------------------------------------------------------------
    ##   x     ... S3 object of class 'Tinflex'
    ##   debug ... if TRUE print intervals
    ## ----------------------------------------------------------------------
    
    genX <- Tinflex.C2R(x)
    Tinflex.print(genX, debug, classname="TinflexC", ...)
}

## --------------------------------------------------------------------------

Tinflex.print <- function(x, debug=FALSE, classname=NULL,...) {
    ## ----------------------------------------------------------------------
    ## S3 method for printing class 'Tinflex' and 'TinflexC'.
    ## ----------------------------------------------------------------------
    
    ## Prepare log-density.
    if (is(x$lpdf, "function")) {
        lpdf.args <- paste(names(formals(x$lpdf)), collapse=",")
        lpdf.body <- gsub("\\s", "", paste(deparse(x$lpdf), collapse=""), perl=TRUE)
        lpdf.body <- gsub("^function\\(\\w+\\)", "", lpdf.body, perl=TRUE)
    } else {
        lpdf.args <- " "
        lpdf.body <- "<none>"
    }
    
    ## Print on console.
    cat("Object of class '",classname,"':\n\n", sep="")
    cat("       log-density(",lpdf.args,") = ",lpdf.body,"\n", sep="")
    cat("       used 2nd derivative: ",x$haved2Tf,"\n\n", sep="")
    
    cat("       area below hat =",x$A.ht,"\n")
    cat("   area below squeeze =",x$A.sq.tot,"\n")
    cat("    ratio hat/squeeze =",x$A.ht/x$A.sq.tot,"   [= upper bound for rejection constant ]\n")
    cat("          # intervals =",ncol(x$ivs)-1,"\n\n")
    
    ## Print boundaries and 'c' values for initial intervals.
    cat(x$iniv,"\n")
    
    if (isTRUE(debug)) {
        ## Print data about for intervals.
        cat("Interval data:\n")
        print(x$ivs)
        cat("\nCumulated areas:\n")
        print(x$Acum)
        cat("\nGuide table:\n")
        print(x$gt)
        cat("\n")
    }

}  ## -- end of Tinflex.print() -- ##

## --------------------------------------------------------------------------
