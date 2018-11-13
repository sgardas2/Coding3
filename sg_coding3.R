myBW = function(x, A, B, w, n.iter = 100){
  # Input: 
  # x: T-by-1 observation sequence
  # A: initial estimate for mz-by-mz transition matrix
  # B: initial estimate for mz-by-mx emission matrix
  # w: initial estimate for mz-by-1 initial distribution over Z_1
  # Output MLE of A and B; we do not update w
  # list(A = A, B=B, w = w)
  
  for(i in 1:n.iter){
    update.para = BW.onestep(x, A, B, w)
    A = update.para$A
    B = update.para$B
    #print(A)
  }
  return(list(A = A, B = B, w = w))
}
BW.onestep = function(x, A, B, w){
  # Input: 
  # x: T-by-1 observation sequence
  # A: current estimate for mz-by-mz transition matrix
  # B: current estimate for mz-by-mx emission matrix
  # w: current estimate for mz-by-1 initial distribution over Z_1
  # Output the updated parameters
  # para = list(A = A1, B = B1)
  
  # We DO NOT update the initial distribution w
  
  T = length(x)
  mz = nrow(A)
  alp = forward.prob(x, A, B, w)
  beta = backward.prob(x, A, B, w)
  myGamma = array(0, dim=c(mz, mz, T-1))
  
  ###
  ## YOUR CODE: 
  ## Compute gamma_t(i,j), which are stored in myGamma
  ##
  for(t in 1:(T-1))
  {
    myGamma[,,t]=t(t(alp[t,]*A)*(B[,x[t+1]]*beta[t+1,]))
  }
  A = rowSums(myGamma, dims = 2)
  A = A/rowSums(A)
  
  tmp = apply(myGamma, c(1, 3), sum)  # mz-by-(T-1)
  tmp = cbind(tmp, colSums(myGamma[, , T-1]))
  for(l in 1:mx){
    B[, l] = rowSums(tmp[, which(x==l)])
  }
  B = B/rowSums(B)
  return(list(A = A, B = B))
}
#You can compute the forward and backward probabilities using the following functions.

forward.prob = function(x, A, B, w){
  
  # Output the forward probability matrix alp 
  # alp: T by mz, (t, i) entry = P(x_{1:t}, Z_t = i)
  
  T = length(x)
  mz = nrow(A)
  alp = matrix(0, T, mz)
  
  # fill in the first row of alp
  alp[1, ] = w*B[, x[1]]
  
  # Recursively compute the remaining rows of alp
  for(t in 2:T){
    tmp = alp[t-1, ] %*% A
    alp[t, ] = tmp * B[, x[t]]
  }
  return(alp)
}

backward.prob = function(x, A, B, w){
  # Output the backward probability matrix beta
  # beta: T by mz, (t, i) entry = P(x_{1:t}, Z_t = i)
  
  T = length(x)
  mz = nrow(A)
  beta = matrix(1, T, mz)
  
  # The last row of beta is all 1.
  # Recursively compute the previous rows of beta
  for(t in (T-1):1){
    tmp = as.matrix(beta[t+1, ] * B[, x[t+1]])  # make tmp a column vector
    beta[t, ] = t(A %*% tmp)
  }
  return(beta)
}


data = read.csv("Coding3_HMM_Data.csv")
dim(data)
mz=2; mx=3
ini.A = matrix(1, mz, mz)
ini.A = ini.A/rowSums(ini.A)
ini.B = matrix(1:6, mz, mx)
ini.B = ini.B/rowSums(ini.B)
ini.w = c(1/2, 1/2)


myout = myBW(data$X, ini.A, ini.B, ini.w, n.iter = 100)


myViterbi=function(x, A, B, w)
{
  d=matrix(0,500,2)
  d[1,]= log(w)+log(B[,x[1]])
  T = length(x)
  for (t in 2:T)
  {
    d[t,]=(apply((d[t-1,])+log(A),2,max))+log(B[,x[t]] )
  }			 
  Z=matrix(,500)
  Z[500]=which.max(d[500,])
  for (t in 499:1)
  {
    Z[t]=which.max(d[t,]+log(A[,Z[t+1]]))
  }
  
  output=matrix(,500)
  for(t in 1:T)
  {
    if(Z[t]==1)
    {k="A"} else {k="B"} 
    output[t]=k
  }
  return(output)
}


##

library(HMM)
hmm0 =initHMM(c("A", "B"), c(1, 2, 3), 
              startProbs = ini.w,
              transProbs = ini.A, emissionProbs = ini.B)

true.out = baumWelch(hmm0, data$X, maxIterations=100, pseudoCount=0)
myout.Z =myViterbi(data$X, true.out$hmm$transProbs, true.out$hmm$emissionProbs, true.out$hmm$startProbs)
true.viterbi = viterbi(true.out$hmm, data$X)
sum(true.viterbi != myout.Z)

myout = myBW(data$X, ini.A, ini.B, ini.w, n.iter = 100)
myout.Z = myViterbi(data$X, myout$A, myout$B, ini.w)
write.table(myout.Z, file = "Coding3_HMM_Viterbi_Output.txt", row.names = FALSE, col.names = FALSE)
A_out=as.data.frame(round(myout$A,7))
colnames(A_out)=c('A','B')
rownames(A_out)=c('A','B')
B_out=as.data.frame(round(myout$B,7))
colnames(B_out)=c('1','2','3')
rownames(B_out)=c('A','B')

pdf("Assignment_Output_3_1061_sgardas2.pdf")
grid.arrange(tableGrob(A_out),tableGrob(B_out),nrow = 2,top = "A and B output after 100 Iterations")

dev.off()


