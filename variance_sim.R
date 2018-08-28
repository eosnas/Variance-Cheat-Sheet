# Simulate variance of population estimates with ratio estimate of observed and estimate detection
varSim <- function(
  M = 1000, #size of strata (number of plots)
  m = 100, #number of sampled plots
  D = 500, #density of items per plot
  mp = 0.5, #detection probability mean
  sp = 0,  #sd of detection probability estimate (se)
  simNum = 10000 #number of simulations
){
  results = matrix(NA, simNum, 6)
  colnames(results) <- c("mean(y)", "sd(y)", "N", "sd1", "sd(mean(y))", "sd(N)")
  p = rnorm(simNum, mp, sp) #detection probability
  N = rpois(M, D)
  for(i in 1:simNum){
    samples = rbinom(m, size=sample(N, m), prob=p[i])
    n = M*mean(samples)/p[i]
    v1 = ((M-m)/M)*var(samples)/m
    v2 = v1 + ((1-p[i])/M)*mean(samples)
    v3 = ((M^2)/(p[i]^2))*(v2 + ((mean(samples)^2)/(p[i]^2)) * var(p))
    results[i,] = c(mean(samples), sd(samples), n, sqrt(v1), sqrt(v2), sqrt(v3))
  }
  return(results)
}

det = seq(0.1, 1, by=0.1)
m = seq(100, 1000, by=100)
simResults = theory1 = theory2 = matrix(NA, length(det), length(m))
for(i in 1:length(m)){for(j in 1:length(det)){
  temp = varSim(m=m[i], mp=det[j])
  simResults[i,j] = sd(temp[,1])
  theory1[i,j] = mean(temp[,4])
  theory2[i,j] = mean(temp[,5])
}}

png("fig1.png")
par(mfrow=c(1,2), pty="s")
plot(m, m, type="n", ylim=c(min(theory2), max(theory2)), xlab="Sample Size", 
     ylab=expression(paste("SD(",bar(y),")")))
for(i in 1:10){
lines(m, theory2[,i], lwd=2, lty=1, col=gray(i/15))
lines(m, simResults[,i], lwd=1, lty=3)
#lines(m, theory1[,i], lwd=1, lty=1, col=gray(i/15))
}
text(x=c(600, 350), y=c(0.15, 1.5), labels = c("p = 0.1", "p = 1.0"))
mtext("A",side=3, adj=0, font=2)
plot(det, det, type="n", ylim=c(0, max(theory2)), xlab="Detection", 
     ylab=expression(paste("SD(",bar(y),")")) )
for(i in 1:10){
  lines(det, theory2[i,], lwd=2, lty=1, col=gray(i/15))
  lines(det, simResults[i,], lwd=1, lty=3)
  #lines(det, theory1[i,], lwd=1, lty=1, col=gray(i/15))
}
text(x=0.5, y=c(0.25, 1.8), labels = c("m = 1000", "m = 100"))
mtext("B",side=3, adj=0, font=2)
par(mfrow=c(1,1), pty="m")
dev.off()

#look at sample vs. binomial variance
#not added to cheat sheet
png("fig2.png")
plot(m, m, type="n", ylim=c(min(theory2), max(theory2)), xlab="Sample Size", 
     ylab=expression(paste("Component of SD(",bar(y),")")))
for(i in 1:10){
  lines(m, theory1[,i], lwd=2, lty=1, col=gray(i/15))
  lines(m, theory2[,i]-theory1[,i], lwd=2, lty=2, col=gray(i/15))
  #lines(m, theory1[,i], lwd=1, lty=1, col=gray(i/15))
}
legend("topright", legend=c("Sample Effort", "Binomial Detection"), lwd=2, lty=c(1,2))
dev.off()
