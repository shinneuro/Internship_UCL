# is it working?

# Probability of Maximum Utility ------------------------------------------

library(OpenMx)
# Matrix A

c0 <- c(1,1,1)
c1 <- c(-1,0,0)
c2 <- c(0,-1,0)
c3 <- c(0,0,-1)

# A1 <- matrix( 
#   c0,c1,c2,c3, nrow=3, ncol=4)
A <- vector("list", length = 4)

A[[1]] <- cbind(c0,c1,c2,c3)
A[[2]] <- cbind(c1,c0,c2,c3)
A[[3]] <- cbind(c1,c2,c0,c3)
A[[4]] <- cbind(c1,c2,c3,c0)

# mean vector
# M_j(t) = A_j*E(t): A[[j]]*E[t,] 

# covariance matrix
# H_j(t) = A_j*diag(S(t) + sg2_epsilon^)*t(A_j)

# P(C(t)=j) = integral(0,inf) multi(Mj(t),H(j(t)))

# vector of probability 

# covariance
for(i in 1:200) {
	for( j in 1:4) {
		cov <- A[[j]]%*%diag(S[i,]+sg2_epsilon)%*%t(A[[j]])
		means <- as.vector(A[[j]] %*% E[i,])
		p[i,j] <- omxMnor(H,M,c(0,0,0),c(Inf,Inf,Inf))
		# or try pmvnorm from the mvtnorm
	}
}