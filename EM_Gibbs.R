####################################################
####################################################
####################################################
# first EM algorithm

len <- c(530, 840, 930) # the lengths
M <- matrix(c(1,1,1,0,1,1,0,0,1), byrow=TRUE, ncol=3) # the transcript sets
k <- c(1666, 896, 81) # the counts for each set
NITER=1024 # no. of iterations
trace_em <- matrix(nrow=ncol(M), ncol=NITER) # save the EM traces here
trace_gibbs <- matrix(nrow=ncol(M), ncol=NITER) # save the Bayesian traces here

#trace_em1=trace_em

# first to find the initial values for EM algorithm
# b is set to the total number of reads (million)
b=sum(k)/1e6
# as the length of each read is 1, thus l=len+1-1. Thus we can use len as l

# from l to calculate s
s=array(0,length(len))
s[1]=len[1]
s[2]=len[2]-s[1]
s[3]=len[3]-s[1]-s[2]

# r.i is a Poisson random varialbe, thus the maximum likelihood estimate of rate is the observable numbers of reads
# initialize
mu.ini=c(k[1]/3,k[1]/3+k[2]/2,k[1]/3+k[2]/2+k[3])/b/len
trace_em[,1]=mu.ini

# EM algorithm implementation
# X is the expected number of reads matrix, i.e., X[i,t] (from the slides page 17)
d=dim(M)
X=array(0,d)
for (iter in 2:NITER){
	print(iter)
	# E step: calcualte the expectation of N
	# i: set index
	for (i in 1:d[1]){
		for (t in 1:d[2]){
		# t: transcript index
			# X[i,t]: nubmer of reads in set i from transcript t
			X[i,t]=k[i]*trace_em[t,iter-1]*M[i,t]/sum(trace_em[,iter-1]*M[i,])
		}
	}
	#print(X)
	# M step: maximization the expectation for E step
	for (index in 1:length(mu)){
		trace_em[index,iter]=sum(X[,index])/sum(b*s*M[,index])
		# another form, where s is not needed:
		#trace_em1[index,iter]=sum(X[,index])/b/len[index]
	}
}

####################################################
####################################################
####################################################
# Bayesian algorithm implementation

library("stats")
# initialize
trace_gibbs[,1]=trace_em[,NITER]
shape=c(1.2,1.2,1.2)
rate=c(0.001,0.001,0.001)
rate=rate+b*len
X=array(0,d)
for (iter in 2:NITER){
print(iter)
# sample from multinomal dist
X[1,]=rmultinom(1,k[1],trace_gibbs[,iter-1]*M[1,]/sum(trace_gibbs[,iter-1]*M[1,]))
X[2,]=rmultinom(1,k[2],trace_gibbs[,iter-1]*M[2,]/sum(trace_gibbs[,iter-1]*M[2,]))
X[3,]=rmultinom(1,k[3],trace_gibbs[,iter-1]*M[3,]/sum(trace_gibbs[,iter-1]*M[3,]))
print(X)
# sample from gamma
shape.new=shape+colSums(X)
trace_gibbs[1,iter]=rgamma(1,shape=shape.new[1],rate=rate[1])
trace_gibbs[2,iter]=rgamma(1,shape=shape.new[2],rate=rate[2])
trace_gibbs[3,iter]=rgamma(1,shape=shape.new[3],rate=rate[3])
}

####################################################
####################################################
####################################################
# plot the result

plot(trace_em[1,],col="red",lty=2,ylim=c(0,1000),xlim=c(0,2000))
lines(trace_em[2,],col="green",lty=2)
lines(trace_em[3,],col="blue",lty=2)

      
lines(trace_gibbs[1,],col="red",lty=1)
lines(trace_gibbs[2,],col="green",lty=1)
lines(trace_gibbs[3,],col="blue",lty=1)

legend("topright",                       
       legend=c("mu1_EM", "mu2_EM", "mu3_EM","mu1_Gibbs","mu2_Gibbs","mu3_Gibbs"),      
       col=c("red", "green", "blue","red", "green", "blue"),    
       lty=c(2,2,2,1,1,1),                    
       lwd=c(1,1,1,1,1,1))


# calcualte the standard error for Gibbs estimates
Gibbs.sd=c(0,0,0)
Gibbs.sd[1]=sd(trace_gibbs[1,]) 
Gibbs.sd[2]=sd(trace_gibbs[2,]) 
Gibbs.sd[3]=sd(trace_gibbs[3,]) 
