Tnw <- function(x, gamma0 = 1, weight = 1, a = 1.0){
	n 	<- length(x)
	if (weight>3) stop("weight must be 1, 2, or 3.")
	Tn 	<- .Call("_Tnw", as.numeric(x), as.integer(n), as.numeric(gamma0), as.numeric(a), as.integer(weight) )
	return(Tn/n)
}
