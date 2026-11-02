# for MS figure 3, adopted from Collatz.R, 10/13/2025

is.even <- function(x) {
  if (x %% 2 == 0) TRUE
  else FALSE
}
is.even(2)
## [1] TRUE
is.even(3)
## [1] FALSE

collatz <- function(n) {
  if (is.even(n)) n/2
  else (3 * n + 1)/2
}

i0<-1
iter <-2000000

out <- matrix(NA, nrow=iter+1-i0, ncol=15)
seq <- matrix(NA, nrow=iter+1-i0, ncol=8)

start_time <- Sys.time()

for (i in i0:iter){

  j<-2*i # for q=111... 
  j<-2*i+1 # for q=110...
  j<-i # for q=11...
  
  n_total <- n <- 4*j-1 # all 12i+3/12i+7/12i+11 = 3(mod4)
  
  q_total <-NULL
  
  n1 <- n
  # while (n != 1)
  while (n>=n1)
  {
    if (is.even(n)) q=0
    else q=1
    q_total <- c(q_total, q)
    n <- collatz(n)
    n_total <- c(n_total, n)
  }
  
  s<-length(q_total)
  r=sum(q_total)
  m = n_total[1]
  Fq<-paste(q_total, collapse='')
  Fq1<-paste('q',Fq, collapse='')
  
  # using geral formula with q, s,r to get F(n) from n
  sum=0
  
  # 1(i) function
  p=gregexpr(pattern ='1',Fq)
  for (k in 1:r)
  {
    sum<-sum+3^(r-k)*2^(p[[1]][k]-1)
  }
  alpha0=r/3
  alpha=alpha0/(1-(2/3)^r)
  Fm=(3^r*m+sum)/2^s
  m1=sum/(2^s-3^r)
  p=r/s
  p2=3^r/2^s
  #Fn
  rem =sum/2^s
  inf =(3^r-2^r)/2^s
  sup =alpha*inf  #observed - updated
  sup1=r*3^(r-1)/2^s #new
  sup2=s*3^r/2^s #from ref [7]
  # sup0=(3^(r-2)*(5+2^(s-r+1))-2^(s-1))/2^s
  
  #Fm/m
  Fmm=n/m
  
  seq[i+1-i0,]<-cbind(s,r,m1,p,p2,m,Fq1,Fmm)
  out[i+1-i0,]<-cbind(i,j,m,n,s,r,m1,p,p2,inf,rem,sup,sup2,Fmm,sup1)
}

end_time <- Sys.time()
end_time-start_time

colnames(seq) <- c('s','r','m1','r/s','3^r/2^s','m','Fq','Fm/m')
colnames(out) <- c('i','j','m','Fm','s','r','m1=sum/(2^s-3^r)','r/s','3^r/2^s','inf','rem','sup','sup2','Fm/m','sup1')

#sort a matrix based on one column
newdata<-out[order(out[,9],decreasing=T),]

Fm <-as.numeric(newdata[,4])
m <-as.numeric(newdata[,3])
s <-as.numeric(newdata[,5])
r <-as.numeric(newdata[,6])
m1 <-as.numeric(newdata[,7])
p2 <-as.numeric(newdata[,9])

# remains with lower and upper bound
inf<-as.numeric(newdata[,10])
rem<-as.numeric(newdata[,11])
sup<-as.numeric(newdata[,12]) #observed - new
sup2<-as.numeric(newdata[,13]) # ref [7]
# sup0<-as.numeric(newdata[,15]) # max
sup1<-as.numeric(newdata[,15]) # New sup

par(mfrow=c(2,2))
#fig4 - A, 600by500  par(mfrow=c(1,1))

plot(p2, rem, frame = T, pch = 20, col = "blue", xlab=expression(paste("Ratio: ",3^r/2^s)), 
     ylab = expression(paste("Second term: ",F[q](m), " and bounds")), 
     type='o',
     # main = "Set: 12i+7",
     main = "A",
     log='y', ylim=c(0.2,800), xlim=c(0.5,1))
lines(p2,sup2,col="black")
lines(p2,inf*max(r)/3,col="red",lty=2)
lines(p2,sup1,col="red")
lines(p2,inf,col="red")

legend("topleft", legend=c("Second term", "Upper bound (Luca, 2005)", "Upper bound (observed)","New upper","Lower bound"),
       col=c("blue", "black","red","red","red"), lty=c(1,1,2,1,1), cex=0.65)


# fig5 - B: ratio Fm/m, 600by500
inf1<-inf/m+p2
rem1<-rem/m+p2
sup0<-sup/m+p2 # new observed
sup3<-sup2/m+p2 #ref [7]
sup4<-sup1/m+p2 # new

a=10
r0=max(r)
a= r0*3^(r0-1)/(3^r0-2^r0)

b=seq(0.5, 1.1, by=0.1)
up<-a+(1-a)*b
low<-1/a+(1-1/a)*b
 
#sup1[sup1>=1]

plot(p2, rem1, frame = T, pch = 20, col = "blue", xlab=expression(paste("Ratio ",3^r/2^s)), 
     ylab = expression(F[q](m)/m), type='o',main = "B", ylim=c(0.5,1.75))
lines(p2,inf1,col="red")
lines(p2,sup4,col="red")
lines(b,up,col="green")
lines(b,low,col="green")
abline(h=1,col="grey")
abline(v=1,col="grey")
text(0.51, 0.62, '1-')
text(0.97, 1.5, '1+')


# fig6 - C: m and m1,557by447
m1_inf <-(3^r-2^r)/(2^s-3^r)
m1_sup <- s*3^r/(2^s-3^r)  # ref [7],   # Ref [7] Florian Luca, On the non-trivial cycles in Collatz’s problem, SUT Journal of Mathematics Vol. 41, No. 1 (2005), 31–41
m1_sup0 <- max(r)/3 * m1_inf  # new observed
m1_sup1 <- r*3^(r-1)/(2^s-3^r)  # new


#jpeg("rplot.jpg", width = 600, height = 500)
plot(p2, m, frame = T, pch = 20, col = "blue", xlab=expression(paste("Ratio ",3^r/2^s)), 
     ylab = expression(paste("Starting number  ",m,", cycle number  ",m[c]," and bounds")),
     type='o', main = "C",log='y', ylim=c(0.15,2100000))
lines(p2,m1, col="green")
lines(p2,m1_sup,col="black")
lines(p2,m1_sup0,col="red",lty=2) #new observed
lines(p2,m1_sup1,col="red")  #new
lines(p2,m1_inf,col="red")
legend("bottomright",legend=c("Starting number", "Cycle number","Upper bound (Luca)",
                          "Upper bound (observed)","New upper","Lower bound"),
       col=c("blue", "green","black","red","red","red"), lty=c(1,1,1,2,1,1), cex=0.65)

 
# 0.83,3.9, 

#fig7  600by500 (fig 3D for laTex files)
# see sr_big.R

options(digits=22)

rm(list = ls())

# install.packages("gmp")
library(gmp)  #for big numbers
library(Rmpfr) #for small numbers

x <- mpfr(306, precBits = 256)/485; print(x)
x <- mpfr(306, precBits = 256)/485; print(x)

precBits <- 256 
one <- mpfr(1, precBits)
two <- mpfr(2, precBits)
three <- mpfr(3, precBits)

# constants used in the computing below 
log32=log(2)/log(3)
log23=log(3)/log(2)

log32e = log(two) / log(three); print(log32e)
log23e = log(three) / log(two); print(log23e)

#method 1 to get numbers (r,s) of sequence 
rmax <-630 #640 this is max of r
# 645 this is the maximum value that we can get in R
alpha=rmax/3 #observed from computing results

out <- matrix(NA, nrow=rmax, ncol=13)
colnames(out) <- c('s','r','r/s','3^r/2^s','inf_m1','a1','log32','Fm/m low',
                   'Fm/m up','m1_max','m1_max2','m1_sup','m1_sup1')

results <- list()

r=1
# s=0
s=(floor(1+r*log23e))

d=as.numeric(s-floor(1+(r-1)/log32e))
q=c(rep(1,r),rep(0,d-1))
q_total <-NULL


start_time <- Sys.time()
start_time 

for (i in 1:rmax){
  r=i
  # s=r
  s=(floor(1+r*log23e)) # new method to find s, 8/15/2025
  
  # this will limit the max of r (645)
  p3r=as.bigz(3)^r
  p2s=as.bigz(2)^s
  p2r=as.bigz(2)^r
  
  p3r1 = as.bigz(3)^(r-1)
  
  b=p3r/p2s
  a=r/s
  a1=(1-1/s)*log32e
  
  m1_inf=((p3r-p2r)/(p2s-p3r))
  m1_sup<-(s*p3r/(p2s-p3r))#ref [Luca, 2005]
  m1_sup1<-(r*p3r1/(p2s-p3r))  #new
  
  # alpha = r/3 / (1-p2r/p3r)
  
  Fm_m_low=1/alpha+(1-1/alpha)*b
  Fm_m_up=alpha+(1-alpha)*b
  
  # 8/15/2025, create sequences w/ max(sum) based on Mike Winkler, 
  # The algorithmic structure of the finite stopping time behavior of 
  # the 3x + 1 function, arXiv:1709.03385 [math.GM], 2017.
  
  d=as.numeric(s-floor(1+(r-1)/log32e))
  ns=as.numeric(s)
  q=c(q,rep(0,d))
  q[ns-d]<-1
  Fq<-paste(q, collapse='')
  Fq1<-paste('q',Fq, collapse='')
  
  # Fq<-paste(q_total, collapse='')
  
  # below from Collatz.R
  # if (is.even(n)) q=0
  # else q=1
  # q_total <- c(q_total, q)
  # Fq<-paste(q_total, collapse='')
  
  # using general formula with q, s,r to get F(n) from n
  sum=0
  
  # 1(i) function
  p=gregexpr(pattern ='1',Fq)
  for (k in 1:r)
  {
    sum<-sum+as.bigz(3)^(r-k)*as.bigz(2)^(p[[1]][k]-1)
  }
  m1_max<-(sum/(p2s-p3r))
  
  # moethd 2: 1(i) = floor(1+(i-1)*log23)
  sum1=0
  for (k in 1:r)
  {
    sum1<-sum1+as.bigz(3)^(r-k)*as.bigz(2)^(floor(1+(k-1)*log23e)-1)
  }
  m1_max2<-(sum1/(p2s-p3r))
  results[[i]] <- list(s=s,r=r,rs=a,p3rp2s=b,m1_inf=m1_inf,a1=a1,log32e=log32e,
                       Fm_m_low=Fm_m_low,Fm_m_up=Fm_m_up,m1_max=m1_max,
                       m1_max2=m1_max2,m1_sup=m1_sup,m1_sup1=m1_sup1)
}

end_time <- Sys.time()
end_time 
end_time-start_time


# Reorder list by x
results <- results[order(sapply(results, function(lst) asNumeric(lst$p3rp2s)))]

# sapply(results, function(lst) format(lst$p3rp2s, digits = 32))

s<-do.call(c, lapply(results, function(z) z$s))
r<-do.call(c, lapply(results, function(z) z$r))
b<-do.call(c, lapply(results, function(z) z$p3rp2s))
m1_sup<-do.call(c, lapply(results, function(z) z$m1_sup))
m1_max<-do.call(c, lapply(results, function(z) z$m1_max))
m1_inf<-do.call(c, lapply(results, function(z) z$m1_inf))
m1_sup1<-do.call(c, lapply(results, function(z) z$m1_sup1))

m1_max2<-do.call(c, lapply(results, function(z) z$m1_max2))

log32e<-do.call(c, lapply(results, function(z) z$log32e))
ymax=max(m1_sup)

# par(mfrow=c(1,1)) 
# ==================
#fig7 - D:  600by500 (fig 3D for laTex files)
#plots of Estimated cycle number bounds
plot(b, m1_sup, frame = T, pch = 20, col = "black", 
     main="D",xlab=expression(paste("Ratio: ",3^r/2^s)), 
     ylab = expression(paste("Estimated bounds of cycle number: ",m[c] )),
     type='l', log='y', ylim=c(0.15,6000000)
)

lines(b,alpha*m1_inf,col="red",lty=2,lwd = 1)
lines(b,m1_sup1,col="red",lty=1)
lines(b,m1_max,col="blue")
lines(b,m1_inf,col="red",lty=1)

abline(v=1,col="black") 

legend("topleft", legend=c("Upper bound (Luca, 2005)",
      paste0("Upper bound (observed, α=",alpha,")"),"new upper bound","Max value", "Lower bound"),
      col=c("black","red","red","blue","red"), lty=c(1,2,1,1), cex=0.8)


 # ======================
# MS Figure 3:
# relationship between m1_sup and m1_max
# conjecture: m1_sup = floor(1+4*log23) * m1_max
# regression:

y<-asNumeric(m1_sup)
x<-asNumeric(m1_max)

# plot(model1)
plot(x, y, frame = T, pch = 20, col = "blue",
     xlab=expression(paste("Calculated Maximum of Cycle number ", m[c])), 
     ylab=expression(paste("Estimated Bound by Luca (2005) for Cycle number ", m[c])), 
     typ='o'
)

model <- lm(y ~ x)
abline(model, col = "red", lwd = 2)

coefs <- coef(model)
eq <- sprintf("y = %.2f + %.2f x", coefs[1], coefs[2])
r2 <- summary(model)$r.squared
r2_txt <- sprintf("R² = %.3f", r2)

# Add text to plot
text(x = min(x), y = max(y), labels = eq, pos = 4, col = "blue")
text(x = min(x), y = 0.95*max(y), labels = r2_txt, pos = 4, col = "blue")

