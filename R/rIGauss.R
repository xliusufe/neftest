rIGauss <- function(n, mu = 1.0, lambda = 1.0){
    rn  = rnorm(n)
    ru  = runif(n)
    x   = .Call("RINV_GAUSS", as.numeric(ru), as.numeric(rn), as.integer(n), as.numeric(mu), as.numeric(lambda))
	return(x)
}
