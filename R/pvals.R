Estimate_shape <- function(x, max.iter, tol){
	ep	= 1e-14
	a	= 0.5
	for(k in 1:max.iter)
	{
		b	= 1/a - trigamma(a)
		c	= log(a) - digamma(a) - log(mean(x)) + mean(log(x))
		if(abs(b)<ep || abs(c/b)<tol)
			break
		else
			a = a - c/b
	}
	a
}

muhat_Gamma <- function(x, alpha){
	beta 	= alpha/mean(x)
	mu		= rep(NA,6)
	for(k in 1:6){
		mu[k] = gamma(k+alpha)/beta^k/gamma(alpha)
	}
	mu
}


pvals <- function(x, distr = "Poisson", bootstrap = FALSE, B = 1000, weight = 1, a = 1.0, max.iter = 100, tol = 1e-8){
	n 	<- length(x)
	if (weight>3) stop("weight must be 1, 2, or 3.")
	if(distr=="Poisson"){
		gamma0 	<- 1
		Tn 		<- .Call("_Tnw", as.numeric(x), as.integer(n), as.numeric(gamma0), as.numeric(a), as.integer(weight) )
		if (bootstrap){
			Tb 			= rep(NA, B)
			lambdahat 	= mean(x)
			for(b in 1:B){
				xb 		<- rpois(n, lambda = lambdahat)
				Tb[b] 	<- .Call("_TnwB", as.numeric(xb), as.integer(n), as.numeric(gamma0), as.numeric(a), as.integer(weight) )
			}
			pval = mean(Tb > Tn)
		}
		else{
			muhat 		= .Call("muhat_poisson", as.numeric(x), as.integer(n))
			sigma2hat 	= .Call("sigma_hat", as.numeric(muhat), as.numeric(gamma0))
			pval 		= pchisq(Tn/sigma2hat, df=1, lower.tail = FALSE)
		}
	}
	else if (distr=="Gamma") {
		gamma0 		<- 2
		Tn 			<- .Call("_Tnw", as.numeric(x), as.integer(n), as.numeric(gamma0), as.numeric(a), as.integer(weight) )
		shapehat 	= Estimate_shape(x, max.iter, tol)
		if (bootstrap){
			Tb 		= rep(NA, B)
			ratehat = shapehat/mean(x)
			for(b in 1:B){
				xb 		<- rgamma(n, shape = shapehat, rate = ratehat)
				Tb[b] 	<- .Call("_TnwB", as.numeric(xb), as.integer(n), as.numeric(gamma0), as.numeric(a), as.integer(weight) )
			}
			pval = mean(Tb > Tn)
		}
		else{
			muhat 		= muhat_Gamma(x, shapehat)
			sigma2hat 	= .Call("sigma_hat", as.numeric(muhat), as.numeric(gamma0))
			pval 		= pchisq(Tn/sigma2hat, df=1, lower.tail = FALSE)
		}
	}
	else if (distr=="Inverse Gaussian") {
		gamma0 	<- 3
		Tn 		<- .Call("_Tnw", as.numeric(x), as.integer(n), as.numeric(gamma0), as.numeric(a), as.integer(weight) )
		if (bootstrap){
			Tb 			= rep(NA, B)
			nuhat 		= mean(x)
			lambdahat 	= 1.0/(mean(1.0/x) - 1.0/nuhat)
			for(b in 1:B){
				xb 		<- rIGauss(n, mu = nuhat, lambda = lambdahat)
				Tb[b] 	<- .Call("_TnwB", as.numeric(xb), as.integer(n), as.numeric(gamma0), as.numeric(a), as.integer(weight) )
			}
			pval = mean(Tb > Tn)
		}
		else{
			muhat 		= .Call("muhat_Ig", as.numeric(x), as.integer(n))
			sigma2hat 	= .Call("sigma_hat", as.numeric(muhat), as.numeric(gamma0))
			pval 		= pchisq(Tn/sigma2hat, df=1, lower.tail = FALSE)
		}
	}
	else
		stop("Distribution must be specified to be one of 'Poisson', 'Gamma', and 'Inverse Gaussian' !")


	return (pval)
}
