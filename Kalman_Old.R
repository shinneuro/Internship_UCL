KalmanPMU_exp <- function(par,data) {
  mLL <- try({
    # parameters
    mu0 <- 0#par[1]
    var0 <- 1000# exp(par[1])
    vart <- exp(par[1])
    vare <- exp(par[2])
    alpha <- exp(par[3])
    lambda <- exp(par[4])
    # setup for fkf
    a0 <- rep(mu0,4)
    P0 <- var0*diag(4)
    dt <- matrix(rep(0,4),ncol=1)
    ct <- matrix(rep(0,4),ncol=1)
    Tt <- array(diag(4),dim=c(4,4,1)) # fixed
    Zt <- array(diag(4),dim=c(4,4,1)) # fixed
    HHt <- array(vart*diag(4),dim=c(4,4,1))
    GGt <- array(vare*diag(4),dim=c(4,4,1))
    yt <- matrix(NA,ncol=4,nrow=nrow(data))
    for(i in 1:4) {
      ut <- data$payoff[data$deck == i]
      ut[ut >= 0] <- ut[ut >= 0]^alpha
      ut[ut < 0] <- -lambda*abs(ut[ut < 0])^alpha
      yt[data$deck == i,i] <- ut
    } 
    filt <- fkf(a0,P0,dt,ct,Tt,Zt,HHt,GGt,t(yt))
    p <- matrix(0,nrow=nrow(data),ncol=4)
    A1 <- matrix(c(1,-1,0,0, 1,0,-1,0, 1,0,0,-1), nrow = 3, byrow = TRUE)
    A <- array(0,dim=c(3,4,4))
    A[,,1] <- A1
    A[,,2] <- A1[,c(2,1,3,4)]
    A[,,3] <- A1[,c(2,3,1,4)]
    A[,,4] <- A1[,c(2,3,4,1)]
    # this part is super slow!
    for(t in 1:nrow(data)) {
      for(i in 1:4) {
        newMu <- as.vector(A[,,i] %*%filt$at[,t])
        newS <- A[,,i] %*% (filt$Pt[,,t] + vare ) %*% t(A[,,i])
        p[t,i] <- pmvnorm(lower=c(0,0,0), mean = newMu, sigma = newS, algorithm=Miwa(steps=128))
        p[p<0] <- 0
      }
    }
    #-2*sum(log(rowSums(cbind(dat$deck==1,dat$deck==2,dat$deck==3,dat$deck==4)*p)))
    p[t(apply(t(filt$at[,1:nrow(data)]),1,function(x) x == max(x)))] <- 0
    rowSums(p)
  })
  if(inherits(mLL,"try-error")) mLL <- NA
  return(mLL)
}