Tnw <- function(x, gamma0 = 1, weight="normal", a = 1.0){
	n 	<- length(x)
	w0 	<-ifelse(weight == "normal", 1, 2)
	Tn 	<- .Call("_Tnw", as.numeric(x), as.integer(n), as.numeric(gamma0), as.numeric(a), as.integer(w0) )
	return(Tn/n)
}
