
hess_ram = function(par,ImpCov,SampCov,Areg,Sreg,A,S,Fmat){

#?sfClusterApplyLB
#sfClusterApplyLB(index, simulationIteration, n, rho, g, h, pv, het=F)

hess.out <- matrix(0,length(par),length(par))
h <- 0.00001


B = solve(diag(nrow(A)) - Areg)
C = diag(nrow(ImpCov)) - solve(ImpCov) %*% SampCov
E = B %*% Sreg %*% t(B)

#if(type=="none"){

# not symmetric

    for(i in 1:nrow(hess.out)){
          for(j in 1:ncol(hess.out)){

          Ai <- (A == i)*1
          Aj <- (A == j)*1
          Aij <- matrix(0,nrow(A),ncol(A))
          Aij[Ai==T & Aj ==T] <- 1

          Si <- (S == i)*1;
          Sj <- (S == j)*1;
          Sij <- matrix(0,nrow(S),ncol(S))
          Sij[Si==T & Sj==T] <- 1

  deriv15_I <- Fmat %*% B %*% Ai %*% E %*% t(Fmat) + Fmat %*% B %*% Si %*% t(B) %*% t(Fmat)
  deriv15_J <- Fmat %*% B %*% Aj %*% E %*% t(Fmat) + Fmat %*% B %*% Sj %*% t(B) %*% t(Fmat)

  # this is cause of asymmetry
  deriv15_IJ <- Fmat %*% B %*% Ai %*% B %*% Aj %*% E %*% t(Fmat) +
                Fmat %*% B %*% Aj %*% B %*% Ai %*% E %*% t(Fmat) +
                Fmat %*% B %*% Ai %*% B %*% Sj %*% t(B) %*% t(Fmat) +
                Fmat %*% B %*% Ai %*% E %*% t(Aj) %*% t(B) %*% t(Fmat) +
                Fmat %*% B %*% Aj %*% B %*% Si %*% t(B) %*% t(Fmat)

  # left out mean part
  hess.out[i,j]  <- trace(solve(ImpCov) %*% deriv15_IJ %*% C - solve(ImpCov) %*%
                          deriv15_I %*% solve(ImpCov) %*% deriv15_J %*% C +
                          solve(ImpCov) %*% deriv15_J %*% solve(ImpCov) %*% deriv15_I %*%
                          solve(ImpCov) %*% SampCov)

          }
      }


hess.out



}
