# get given r,s and q to get m


source("C:/Users/Owner/OneDrive/3n+1/Collatz_environment.R")


# generate a string with max summary:
# r=19, s=31, 1101101101011011010110110110100


maxq_list <- list()
q <-'10' # Initial of q 
for (r in 2:20){
  s=floor(1+r*log23e)
  s1=floor(1+(r-1)*log23e)
  d=as.numeric(s-s1)
  q<-paste0(substr(q, 1, nchar(q) - 1), "1") #replace last digit '0' with '1'
  zeros=paste0(rep(0,d), collapse='') # generate numbers of '0's with length of d  
  q=paste0(q, zeros) # added numbers of '0's with length of d to the q
  maxq_list[[r-1]]<-list(r,s,q)
}


# find m with max summation from given q
start_time <- Sys.time()
maxq_list <- list()
for (r in 1000){
  q <-'10' # Initial of q 
  for (r1 in 2:r){
    # generate a string with max summary: e.g.
    # r=19, s=31, 1101101101011011010110110110100
    s=as.numeric(floor(1+r1*log23e))
    s1=as.numeric(floor(1+(r1-1)*log23e))
    d=(s-s1)
  
    #replace last digit '0' with '1'
    q<-paste0(substr(q, 1, nchar(q) - 1), "1")
    # generate numbers of '0's with length of d 
    zeros=paste0(rep(0,d), collapse='')
    # added numbers of '0's with length of d to the q
    q=paste0(q, zeros)
  }
  print(c(r,s,q))

  # initialize the values and integer type (big z) 
  
  m=-as.bigz(1)
  Fm=-as.bigz(1)
  p=as.bigz(2)
  t=as.bigz(2)

  for (i in 3:s){
    step <- substr(q, i, i)
    if (step==0) {
      p=p #power of 3, max=r
      if (is.even(Fm)) {
        Fm=Fm/2
        m=m
      } else {
        Fm=(3^p+Fm)/2
        m=2^t+m
      }
    } else if (step==1){
      p=p+1
      if (is.even(Fm)) {
        Fm=(3*(3^(p-1)+Fm)+1)/2
        m=2^t+m
      } else {
        Fm=(3*Fm+1)/2
        m=m
      }
    }
    t=t+1 # power of 2 max=s   
  }
  maxq_list[[r-1]]<-list(r=r,s=s,q=q,m=m,Fm=Fm)
}
end_time <- Sys.time()
end_time-start_time

m
Fm
t
p
q

# get all m's
M=NULL
FM=NULL
Q=NULL
R=NULL
S=NULL
for (i in 1:length(maxq_list)){
  M[i]=maxq_list[[i]]$m
  FM[i]=maxq_list[[i]]$Fm
  Q[i]=maxq_list[[i]]$q
  R[i]=maxq_list[[i]]$r
  S[i]=maxq_list[[i]]$s
}



# min summation q
#q='1111111111111111111000000000000'

minq_list <- list()
start_time <- Sys.time()
for (r in 1000) {

# s=floor(1+r*log2(3))
s=as.numeric(floor(1+r*log23e))

ones=paste0(rep(1,r), collapse='')  
zeros=paste0(rep(0,s-r), collapse='') 
q=paste0(ones, zeros)

# m=-one
# Fm=-one
m=-as.bigz(1)
Fm=-as.bigz(1)
p=as.bigz(2)
t=as.bigz(2)

 # initialize the values and integer type (big z) 
  
  m=-as.bigz(1)
  Fm=-as.bigz(1)
  p=as.bigz(2)
  t=as.bigz(2)

  for (i in 3:s){
    step <- substr(q, i, i)
    if (step==0) {
      p=p #power of 3, max=r
      if (is.even(Fm)) {
        Fm=Fm/2
        m=m
      } else {
        Fm=(3^p+Fm)/2
        m=2^t+m
      }
    } else if (step==1){
      p=p+1
      if (is.even(Fm)) {
        Fm=(3*(3^(p-1)+Fm)+1)/2
        m=2^t+m
      } else {
        Fm=(3*Fm+1)/2
        m=m
      }
    }
    t=t+1 # power of 2 max=s   
  }
# odd number
odd=(m+1)/(as.bigz(2)^r)
minq_list[[r-1]]<-list(r=r,s=s,q=q,m=m,Fm=Fm,odd=odd)
}
end_time <- Sys.time()
end_time-start_time



# get all m's
M=NULL
FM=NULL
Q=NULL
R=NULL
S=NULL
ODD=NULL
for (i in 1:length(minq_list)){
  M[i]=minq_list[[i]]$m
  FM[i]=minq_list[[i]]$Fm
  Q[i]=minq_list[[i]]$q
  R[i]=minq_list[[i]]$r
  S[i]=minq_list[[i]]$s
  ODD[i]=minq_list[[i]]$odd
}

