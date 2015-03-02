## Given the components of MME, estimate variance components by REML in the Newton Raphson framework.

## Arguments
## A: an additive relationship matrix
## y: a vector of response
## X: a incidence matrix for fixed effects
## Z: a incidence matrix for random effects
## initE: initial value for the residual variance
## initU: initial value for the additive genetic variance
## disp: a logical value. If true, display estimates of variance components at each iteration.

`nrReml` <-
function(A, y, X, Z, initE, initU, disp){

	N <- length(y)
	Ze <- diag(N)
	Var <- c(initU, initE)
	H <- matrix(0, ncol=2, nrow=2)
	s <- matrix(0, ncol=1, nrow=2)
	diff1 <- 10
	diff2 <- 10

	i <- 0
	while (diff1 > 10e-7 & diff2 > 10e-7 ){
	
		i <- i + 1
	
		G <- A*Var[1]
		R <- Ze%*%t(Ze)*Var[2]
		V <- Z%*%G%*%t(Z) + R
		Vinv <- solve(V)
		P <- Vinv - Vinv%*%X%*%solve(t(X)%*%Vinv%*%X)%*%t(X)%*%Vinv

		H[1,1] <- -sum(diag((P%*%Z%*%t(Z)%*%P%*%Z%*%t(Z) ))) + 2*t(y)%*%P%*%Z%*%t(Z)%*%P%*%Z%*%t(Z)%*%P%*%y
		H[1,2] <- -sum(diag((P%*%Z%*%t(Z)%*%P%*%Ze%*%t(Ze) ))) + 2*t(y)%*%P%*%Z%*%t(Z)%*%P%*%Ze%*%t(Ze)%*%P%*%y

		H[2,1] <- H[1,2]
		H[2,2] <- -sum(diag((P%*%Ze%*%t(Ze)%*%P%*%Ze%*%t(Ze) ))) + 2*t(y)%*%P%*%Ze%*%t(Ze)%*%P%*%Ze%*%t(Ze)%*%P%*%y
	
		s[1,1] <- sum(diag((P%*%Z%*%t(Z) ))) - (t(y)%*%P%*%Z%*%t(Z)%*%P%*%y )
		s[2,1] <- sum(diag((P%*%Ze%*%t(Ze) ))) - (t(y)%*%P%*%Ze%*%t(Ze)%*%P%*%y )
	 
		newVar <- Var - solve(H)%*%s
		diff1 <- abs(Var[1] - newVar[1])
		diff2 <- abs(Var[2] - newVar[2])
		Var <- newVar
		if (disp == TRUE){
			cat('\n')
			cat("iteration ", i, '\n')
			cat("sig2U", Var[1], '\n')		
			cat("sig2E", Var[2], '\n')
		}
	
	}
	
	cat('\n')
	return(list( "sigma2U"=Var[1], "sigma2E"=Var[2]   ))

}




## Given the components of MME, estimate variance components by REML in the Fisher's Scoring framework.

## Arguments
## A: an additive relationship matrix
## y: a vector of response
## X: a incidence matrix for fixed effects
## Z: a incidence matrix for random effects
## initE: initial value for the residual variance
## initU: initial value for the additive genetic variance
## disp: a logical value. If true, display estimates of variance components at each iteration.

`fsReml` <-
function(A, y, X, Z, initE, initU, disp){

	N <- length(y)
	Ze <- diag(N)
	Var <- c(initU, initE)
	I <- matrix(0, ncol=2, nrow=2)
	s <- matrix(0, ncol=1, nrow=2)
	diff1 <- 10
	diff2 <- 10

	i <- 0
	while (diff1 > 10e-7 & diff2 > 10e-7 ){
	
		i <- i + 1
	
		G <- A*Var[1]
		R <- Ze%*%t(Ze)*Var[2]
		V <- Z%*%G%*%t(Z) + R
		Vinv <- solve(V)
		P <- Vinv - Vinv%*%X%*%solve(t(X)%*%Vinv%*%X)%*%t(X)%*%Vinv

		I[1,1] <- sum(diag((P%*%Z%*%t(Z)%*%P%*%Z%*%t(Z) )))
		I[1,2] <- sum(diag((P%*%Z%*%t(Z)%*%P%*%Ze%*%t(Ze) ))) 
		I[2,1] <- H[1,2]
		I[2,2] <- sum(diag((P%*%Ze%*%t(Ze)%*%P%*%Ze%*%t(Ze) ))) 	
		s[1,1] <- sum(diag((P%*%Z%*%t(Z) ))) - (t(y)%*%P%*%Z%*%t(Z)%*%P%*%y )
		s[2,1] <- sum(diag((P%*%Ze%*%t(Ze) ))) - (t(y)%*%P%*%Ze%*%t(Ze)%*%P%*%y )
	 
		newVar <- Var - solve(H)%*%s
		diff1 <- abs(Var[1] - newVar[1])
		diff2 <- abs(Var[2] - newVar[2])
		Var <- newVar
		if (disp == TRUE){
			cat('\n')
			cat("iteration ", i, '\n')
			cat("sig2U", Var[1], '\n')		
			cat("sig2E", Var[2], '\n')
		}
	
	}
	
	cat('\n')
	return(list( "sigma2U"=Var[1], "sigma2E"=Var[2]   ))

}




## Given the components of MME, estimate variance components by EM Algorithm. 

## Arguments
## Ainv: an inverse of additive relationship matrix
## y: a vector of response
## X: a incidence matrix for fixed effects
## Z: a incidence matrix for random effects
## initE: initial value for the residual variance
## initU: initial value for the additive genetic variance
## disp: a logical value. If true, display estimates of variance components at each iteration.

`emreml` <-
function(Ainv, y, X, Z, initE, initU, disp){

	# MME
	n <- length(y)
	Xpy <- t(X)%*%y
	Zpy <- t(Z)%*%y
	XpX <- t(X)%*% X
	XpZ <- t(X)%*%Z
	ZpX <- t(Z)%*%X
	ZpZ <- t(Z)%*%Z
	RHS <- c(Xpy, Zpy)
	oldE <- initE
	oldU <- initU
	rankX <- qr(X,LAPACK=TRUE)$rank
	rankA <- qr(Ainv,LAPACK=TRUE)$rank
	
	lhsRow <- length(RHS)
	z <- length(RHS) - dim(ZpZ)[1] + 1

	diff1 <- 1
	diff2 <- 1
	i <- 0
	while (diff1 > 10E-6 & diff2 > 10E-6){
		i <- i+1
		alpha <- as.vector((oldE/oldU))
		LHS <- rbind( cbind(XpX, XpZ), cbind(ZpX, ZpZ+Ainv*alpha ) )
		B <- solve(LHS)%*%RHS
		e <-  y - cbind(X,Z)%*%B
		sig2E <- (t(e)%*%y)/(n-rankX)
		c22 <- solve(LHS)[z:lhsRow, z:lhsRow]
		u <- B[z:length(B)]
		# sum(Ainv*c22) is same as sum(diag(Ainv%*%c22))
		sig2U <- (t(u)%*%Ainv%*%u + sum(Ainv*c22)*oldE )/(rankA)
		diff1 <- abs(sig2E-oldE)
		diff2 <- abs(sig2U-oldU)
		if (disp == TRUE){
			cat('\n')
			cat("iteration ", i, '\n')
			#cat("oldE", oldE, '\n')
			cat("sig2E", sig2E, '\n')
			#cat("oldU", oldU, '\n')
			cat("sig2U", sig2U, '\n')				
		}
		oldE <- sig2E
		oldU <- sig2U
	}
	
	cat('\n')
	return(list( "sigma2E"=sig2E, "sigma2U"=sig2U   ))

}





## Given the components of MME, estimate variance components by REML in the average information framework.

## Arguments
## A: an additive relationship matrix
## y: a vector of response
## X: a incidence matrix for fixed effects
## Z: a incidence matrix for random effects
## initE: initial value for the residual variance
## initU: initial value for the additive genetic variance
## disp: a logical value. If true, display estimates of variance components at each iteration.

`aiReml` <-
function(A, y, X, Z, initE, initU, disp){

	N <- length(y)
	Ze <- diag(N)
	Var <- c(initU, initE)
	AI <- matrix(0, ncol=2, nrow=2)
	s <- matrix(0, ncol=1, nrow=2)
	diff1 <- 10
	diff2 <- 10

	i <- 0
	while (diff1 > 10e-7 & diff2 > 10e-7 ){
	
		i <- i + 1
	
		G <- A*Var[1]
		R <- Ze%*%t(Ze)*Var[2]
		V <- Z%*%G%*%t(Z) + R
		Vinv <- solve(V)
		P <- Vinv - Vinv%*%X%*%solve(t(X)%*%Vinv%*%X)%*%t(X)%*%Vinv

		AI[1,1] <- sum(diag((t(y)%*%P%*%Z%*%t(Z)%*%P%*%Z%*%t(Z)%*%P%*%y )))
		AI[1,2] <-  sum(diag((t(y)%*%P%*%Z%*%t(Z)%*%P%*%Ze%*%t(Ze)%*%P%*%y ))) 
		AI[2,1] <- H[1,2]
		AI[2,2] <- sum(diag((t(y)%*%P%*%Ze%*%t(Ze)%*%P%*%Ze%*%t(Ze)%*%P%*%y ))) 	
		s[1,1] <- sum(diag((P%*%Z%*%t(Z) ))) - (t(y)%*%P%*%Z%*%t(Z)%*%P%*%y )
		s[2,1] <- sum(diag((P%*%Ze%*%t(Ze) ))) - (t(y)%*%P%*%Ze%*%t(Ze)%*%P%*%y )
	 
		newVar <- Var - solve(AI)%*%s
		diff1 <- abs(Var[1] - newVar[1])
		diff2 <- abs(Var[2] - newVar[2])
		Var <- newVar
		if (disp == TRUE){
			cat('\n')
			cat("iteration ", i, '\n')
			cat("sig2U", Var[1], '\n')		
			cat("sig2E", Var[2], '\n')
		}
	
	}
	
	cat('\n')
	return(list( "sigma2U"=Var[1], "sigma2E"=Var[2]   ))

}
