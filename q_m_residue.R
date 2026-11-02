# given r,s to get q, m and no. of residue (https://oeis.org/A100982): m(mod2^s) in Collartz process

# based on the sequence tree:

# no. of residue (https://oeis.org/A100982): m(mod2^s) = sum(number of 0's in the tails of previous q's)

 
# how grow a tree branch:
# add 10 before each "0" if d=2 else,
# add 10 before each "0" if d=1.
                         
# e.g.,     d=1       d=2
# start 1100 -> 11010 -> 1101100
#             | 11100 -> 1110100
#                        1111000

# 10/25/2025

source("C:/Users/Owner/OneDrive/3n+1/Collatz_environment.R")

# function to insert a substring in a string at certain position
insert_substr <- function(orig_str, substr, pos) {
  # Split the original string into two parts
  part1 <- substr(orig_str, 1, pos - 1)
  part2 <- substr(orig_str, pos, nchar(orig_str))
  
  # Concatenate the parts with the substring in between
  result <- paste0(part1, substr, part2)
  return(result)
}



start_time <- Sys.time()
start_time

# initialize                     
e<-NULL
n<- NULL
E<- NULL
r=2
s=floor(1+r*log23)
q <- "1100"
k=1
e[k]=2  # # of zeros of the sequence tail 
E[[r]]$e<-e
E[[r]]$q<-q

# not necessary
E[[r]]$m<-3
E[[r]]$d<-2
E[[r]]$r<-2
E[[r]]$s<-4

E[[1]]$r<-1
E[[1]]$s<-2
E[[1]]$e<-1
E[[1]]$q<-'10'
E[[1]]$m<-1
E[[1]]$d<-1

mk<-NULL
Fmk<-NULL

# start from 3:29 took 3.698423559202088206632 hours
for (r in 3:12){
  s=floor(1+r*log23)
  s1=floor(1+(r-1)*log23)
  d=s-s1
  ins <- ifelse(d == 2, "10", "1")
  
  pre_q <- E[[r-1]]$q
  e1<-E[[r-1]]$e

  n <- sum(e) # No. of residue classes in m(mod2^s)
  
  k=1
  # for each previous q's'
  for (u in 1:length(e1)){
    j=d  
    q0<-pre_q[u]
    # get new q's  
    for (i in 1:e1[u]){
      position <- s1-i+1
      q[k] <- insert_substr(q0, ins, position)  
      # cat("new q",k,"=",q[k], "\n")
      e[k] <- j
      j=j+1

      # get m and Fm from q created above
      # initialize the values and integer type (big z) 
    
      m=-1
      Fm=-1
      p=2
      t=2
      
      # use when s>50, about 3X slow
      # m=-big_1
      # Fm=-big_1
      # p=big_2
      # t=big_2
    
      for (i in 3:s){
        step <- substr(q[k], i, i)
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
    mk[k]<-m
    Fmk[k]<-Fm
    
    # Using big number package  
    # mk[k]<-as.numeric(m)
    # Fmk[k]<-as.numeric(Fm)
    k=k+1
    }
  }
  E[[r]]<-list(r=r,s=s,d=d,e=e,n=n,q=q, m=mk,Fm=Fmk)
}
end_time <- Sys.time()
start_time
end_time
end_time-start_time
#End of temp




# get all m's
M=NULL
FM=NULL
Q=NULL
R=NULL
S=NULL
N=NULL
D=NULL

for (i in 2:length(E)){
  M[[i]]=E[[i]]$m
  FM[[i]]=E[[i]]$Fm
  Q[[i]]=E[[i]]$q
  R[i]=E[[i]]$r
  S[i]=E[[i]]$s
  N[i]=E[[i]]$n
  D[i]=E[[i]]$d
}
