
# the second column is the number of residual classes of mod(2^s), r<=13, s<=21
# 8/13/2025

file_path <- "https://oeis.org/A177789/b177789.txt"

data <- read.table(file_path, header = FALSE)
head(data)

# Remove the first two rows
data <- data[-c(1, 2), ]

# Rename the column
colnames(data) <- c('j','m')

# search range for
# constants:
log23=log(3)/log(2)
log32=log(2)/log(3)


r0=20
s0=floor(1+r0*log23);

i0<-0
iter <-length(data[,1])-1
# iter <-floor((2^s0-3)/4)
# 4*iter+3 = 2^s0

out <- matrix(NA, nrow=iter+1-i0, ncol=15)
seq <- matrix(NA, nrow=iter+1-i0, ncol=9)

start_time <- Sys.time()
for (i in i0:iter){
  # j<-2*(i)+1
  # n_total <- n <- 6*j-3  #12i+3
  # n_total <- n <- 6*j+1  #12i+7
  # n_total <- n <- 6*j+5  #12i+11
  
  j<-i
  # n_total <- n <- 4*j+3 # all 12i+3/12i+7/12i+11
  n_total <- n <- data[j+1,2]
  
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
  
  alpha=10
  p3r=3^r
  p2s=2^s
  p2r=2^r
  
  Fm=(3^r*m+sum)/p2s
  m1=sum/(p2s-p3r)
  p=r/s
  p2=p3r/p2s
  #Fn
  rem =sum/2^s
  inf =(p3r-p2r)/p2s
  sup =alpha*inf  #observed
  sup2=s*p3r/p2s #from ref [7]
  
  # sup0=(3^(r-2)*(5+2^(s-r+1))-2^(s-1))/p2s
  
  # moethd 2: 1(i) = floor(1+(i-1)*log23)
  sum1=0
  for (k in 1:r)
  {
    sum1<-sum1+3^(r-k)*2^(floor(1+(k-1)*log23)-1)
  }
  m1_max<-sum1/(p2s-p3r)
  
  
  #Fm/m
  Fmm=n/m
  
  seq[i+1-i0,]<-cbind(s,r,m1,p,p2,m,Fq1,Fmm,m1_max)
  out[i+1-i0,]<-cbind(i,j,m,n,s,r,m1,p,p2,inf,rem,sup,sup2,Fmm,m1_max)
}

end_time <- Sys.time()
end_time-start_time

colnames(seq) <- c('s','r','m1','r/s','3^r/2^s','m','Fq','Fm/m','m1_max')
colnames(out) <- c('i','j','m','Fm','s','r','m1=sum/(2^s-3^r)','r/s','3^r/2^s','inf','rem','sup','sup2','Fm/m','m1_max')
#sort a matrix based on one column
newdata<-out[order(out[,9],decreasing=T),]

Fm <-as.numeric(newdata[,4])
m <-as.numeric(newdata[,3])
s <-as.numeric(newdata[,5])
r <-as.numeric(newdata[,6])
m1 <-as.numeric(newdata[,7])
p2 <-as.numeric(newdata[,9])

m1_inf <-(p3r-p2r)/(p2s-p3r)
m1_sup <- s*p3r/(p2s-p3r)  # ref [7]


# remains with lower and upper bound
inf<-as.numeric(newdata[,10])
rem<-as.numeric(newdata[,11])
sup<-as.numeric(newdata[,12]) #observed
sup2<-as.numeric(newdata[,13]) # ref [7]
m1_max<-as.numeric(newdata[,15]) # max

# find 11...10...0 sequence number of 1 = r0
q10<-(seq[seq[,2]==13,])
q11<-q10[order(q10[,7],decreasing=F),]



# check max of k-cycle / m-cycle: in {He23}
r=306
r=431166034846567
r=160
s=floor(1+r*log23)
k=floor(1+(r-1)*log23)-r+1




