options(digits=22)

rm(list = ls())

# constants used in the computing below 
log32=log(2)/log(3)
log23=log(3)/log(2)

#method 1 to get numbers (r,s) of sequence 
rmax <-640 #this is max of r, 645 this is the maximum value that we can get in R
alpha=rmax/4 #observed from computing results

out <- matrix(NA, nrow=rmax, ncol=12)
colnames(out) <- c('s','r','r/s','3^r/2^s','inf_m1','a1','log32','Fm/m low','Fm/m up','m1_max','m1_max2','m1_sup')

r=1
# s=0
s=floor(1+r/log32)

d=s-floor(1+(r-1)/log32)
q=c(rep(1,r),rep(0,d))
q_total <-NULL

for (i in 2:rmax){
r=i
# s=r
s=floor(1+r/log32) # new method to find s, 8/15/2025

# this will limit the max of r (645)
p3r=3^r
p2s=2^s
p2r=2^r

b=p3r/p2s

# old method to find s by using 1/2 < 3^r/2^s < 1
# while (b<=0.5 | b>=1)
# {
# b=3^r/2^s
# s<-s+1
# }
# s<-s-1

a=r/s

a1=(1-1/s)*log32

m1_inf=(p3r-p2r)/(p2s-p3r)
m1_sup<-s*p3r/(p2s-p3r) #ref [Luca, 2005]

Fm_m_low=1/alpha+(1-1/alpha)*b
Fm_m_up=alpha+(1-alpha)*b

# 8/15/2025, create sequences w/ max(sum) based on Mike Winkler, The algorithmic structure of the finite stopping time behavior of the 3x + 1 function, arXiv:1709.03385 [math.GM], 2017.

d=s-floor(1+(r-1)/log32)
q=c(q,rep(0,d))
q[s-d]<-1
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
  sum<-sum+3^(r-k)*2^(p[[1]][k]-1)
}
m1_max<-sum/(p2s-p3r)

# moethd 2: 1(i) = floor(1+(i-1)*log23)
sum1=0
for (k in 1:r)
{
  sum1<-sum1+3^(r-k)*2^(floor(1+(k-1)*log23)-1)
}
m1_max2<-sum1/(p2s-p3r)

out[i,]<-cbind(s,r,a,b,m1_inf,a1,log32,Fm_m_low,Fm_m_up,m1_max,m1_max2,m1_sup)
}

out1<-out[order(out[,4],decreasing=F),]
# remove any row with NA/inf
sortout <- out1[apply(out1, 1, function(row) all(is.finite(row))), , drop=FALSE]

r <-as.numeric(sortout[,2])
s <-as.numeric(sortout[,1])
a <-as.numeric(sortout[,3])
b <-as.numeric(sortout[,4])
m1_inf <-as.numeric(sortout[,5])
a1 <-as.numeric(sortout[,6])
a2 <-as.numeric(sortout[,7])


Fm_m_low<-as.numeric(sortout[,8])
Fm_m_up<-as.numeric(sortout[,9])

m1_max<-as.numeric(sortout[,10])
m1_max2<-as.numeric(sortout[,11])
m1_sup<-as.numeric(sortout[,12])

max(m1_sup)

# par(mfrow=c(1,1)) 
# ==================
#fig7  600by500 (fig 3D for laTex files)
#plots of Estimated cycle number bounds
plot(b, m1_sup, frame = T, pch = 20, col = "black", 
     main="D",xlab=expression(paste("Ratio ",3^r/2^s)), 
     ylab = expression(paste("Estimated bounds of cycle number: ",m[c] )),
     type='l', log='y', ylim=c(0.15,1000000)
     #ylab = expression(paste("Estimated cycle number bounds")), type='o', log='y'),
)

lines(b,alpha*m1_inf,col="red",lty=2,lwd = 1)
lines(b,m1_max,col="blue")
lines(b,m1_inf,col="red",lty=1)

abline(v=1,col="black")

legend("topleft", legend=c("Upper bound (Luca, 2005)", paste0("Upper bound (observed, α=",alpha,")"),"Max value", "Lower bound"),
       col=c("black", "red","blue","red"), lty=c(1,2,1,1), cex=0.8)

# ======================
# relationship between m1_sup and m1_max
# conjecture: m1_sup = floor(1+4*log23) * m1_max
# regression:

y<-(m1_sup)
x<-(m1_max)

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

# =============

y<-r
x<-b


# plot(model1)
plot(x, y, frame = T, pch = 20, col = "blue",
     xlab=expression(paste("Calculated Maximum of Cycle number ", m[c])), 
     ylab=expression(paste("Estimated Bound by Luca (2005) for Cycle number ", m[c])), 
     typ='o'
)


y1<-r[seq(638, 638-11, by = -2)]
x1<-b[seq(638, 638-11, by = -2)]

model <- lm(y1 ~ x1)
abline(model, col = "red", lwd = 2)
coefs <- coef(model)

for (i in 0:10){
  y1<-r[seq(637-i*12, 637-11-i*12, by = -2)]
  x1<-b[seq(637-i*12, 637-11-i*12, by = -2)]
  # print(y1)
  model <- lm(y1 ~ x1)
  abline(model, col = "red", lwd = 2)
  print(floor(sum(coef(model))))
  
}


y1<-r[seq(1, 12, by = 2)]
x1<-b[seq(1, 12, by = 2)]

y1<-r[seq(2, 24, by = 2)]
x1<-b[seq(2, 24, by = 2)]

y1<-r[seq(13, 36, by = 2)]
x1<-b[seq(13, 36, by = 2)]

y1<-c(r[seq(26, 36, by = 2)], r[seq(37, 48, by = 2)])
x1<-c(b[seq(26, 36, by = 2)], b[seq(37, 48, by = 2)])

model <- lm(y1 ~ x1)
abline(model, col = "red", lwd = 2)
coefs <- coef(model)
# predict r when 3^r/2^s=1
print(floor(sum(coef(model))))
# 


# ==========







which.max(m1_inf)

which.max(a)
which.max(b)

# 6/11/2025
#new plot of Fm/m with boudaries 
plot(b, Fm_m_low, frame = T, pch = 20, col = "black", 
     main="D",xlab=expression(paste("Ratio ",3^r/2^s)), 
     ylab = expression(paste("Boundaries of Fm/m")),
     type='l', ylim=c(0.5,6)
)
lines(b,Fm_m_up,col="blue",lty=c(1))
abline(h = 1, col = "red")



#4/28/2025, alternative fig2

sortout<-out[order(out[,3],decreasing=F),]
r <-as.numeric(sortout[,2])
s <-as.numeric(sortout[,1])
a <-as.numeric(sortout[,3])
b <-as.numeric(sortout[,4])
m1_inf <-as.numeric(sortout[,5])
a1 <-as.numeric(sortout[,6])
a2 <-as.numeric(sortout[,7])
m1_sup<-s*3^r/(2^s-3^r) #ref [7]

max(m1_sup)


#plots of Estimated ratio of r and s
plot(b, a, frame = T, pch = 20, col = "blue",  xlab=expression(paste("Ratio ",3^r/2^s)), 
     ylab =expression(paste("Ratio of  ",r/s)), type='o'
     #    ,xlim=c(0.85,1)
     #   ,ylim=c(0.625,0.63095)
)

plot(a, b, frame = T, pch = 20, col = "blue", type='o'
     ,ylab=expression(paste("Ratio ",3^r/2^s))
     ,xlab =expression(paste("Ratio of  ",r/s))
     ,ylim=c(0.39,1)
     ,xlim=c(0.6230,0.63095)
)

par(mfrow=c(1,2)) 
colors <- rainbow(length(a)) 
plot(a, b, frame = F, pch = 20,  col = colors, xlab=expression(paste("Ratio of  ",r/s)), 
     ylab =expression(paste("Ratio of ",3^r/2^s)), type='p',xlim=c(0.6305,0.631),ylim=c(0.4,1), main="A")
lines(a, b,col="black") 

plot(a, b, frame = F, pch = 20,  col = colors, xlab=expression(paste("Ratio of  ",r/s)), 
     ylab =expression(paste("Ratio of ",3^r/2^s)), type='p',xlim=c(0.63,0.631),ylim=c(0.4,1), main="B")




#plot of "Ratio of ",3^r/2^s vs minimum steps (s) 

sortout<-out[order(out[,1],decreasing=T),]
r <-as.numeric(sortout[,2])
s <-as.numeric(sortout[,1])
a <-as.numeric(sortout[,3])
b <-as.numeric(sortout[,4])
m1_inf <-as.numeric(sortout[,5])
a1 <-as.numeric(sortout[,6])
a2 <-as.numeric(sortout[,7])

# will combine fig 2,3 and save as fig2.PNG
par(mfrow=c(1,2))  

#fig2
plot(s, a, frame = T, pch = 20, col = "blue", xlab=expression(paste("minimum steps (s)")), 
     ylab =expression(paste("Ratio of  ",r/s)), type='o',xlim=c(0,100),ylim=c(0.3,0.65), main="A")
lines(s,a1,col="red")
lines(s,a2,col="red")


plot(s, b, frame = T, pch = 20, col = "blue", xlab=expression(paste("minimum steps (s)")), 
     ylab =expression(paste("Ratio of ",3^r/2^s)), type='o',xlim=c(0,100),ylim=c(0.3,1), main="B")
lines(s,a/a,col="red")
lines(s,0.5*a/a,col="red")

#4/28/2025

par(mfrow=c(1,2)) 
colors <- rainbow(length(a)) 
plot(a, b, frame = F, pch = 20,  col = colors, xlab=expression(paste("Ratio of  ",r/s)), 
     ylab =expression(paste("Ratio of ",3^r/2^s)), type='p',xlim=c(0.625,0.631),ylim=c(0.4,1), main="A")

plot(a, b, frame = F, pch = 20,  col = colors, xlab=expression(paste("Ratio of  ",r/s)), 
     ylab =expression(paste("Ratio of ",3^r/2^s)), type='p',xlim=c(0.63,0.631),ylim=c(0.4,1), main="B")




#6/24/2025
#method 2  to get numbers (r,s) of uence, get Table 4 numbers
#	Mike Winkler, On the structure and the behaviour of Collatz 3n + 1 sequences, arXiv:1412.0519v1 [math.GM] 27 Nov 2014
# avoid calculating 3^r/2^s
# using formula: s=floor(1+r/a2)
start_time <- Sys.time()
#r range: 
r1=1
r2=2000
out <- matrix(NA, nrow=r2-r1+1, ncol=5)
colnames(out) <- c('s','r','r/s','a1','a2')

r=0
s=0
a2=log(2)/log(3)

for (i in r1:r2){
  r=i
  s=floor(1+r/a2)
  a=r/s
  a1=(1-1/s)*a2
  
  #  m1_inf=(3^(r-1)-2^(r-1))/(2^s-3^r)
  out[i-r1+1,]<-cbind(s,r,a,a1,a2)
}

end_time <- Sys.time()
end_time-start_time

#c=0.630927835
#c=0.63092975
c=0

out1<-out[which(out[,3]>=c),]
sortout<-out1[order(out1[,1],decreasing=F),]

s <-as.numeric(sortout[,1])
r <-as.numeric(sortout[,2])
a <-as.numeric(sortout[,3])
a1 <-as.numeric(sortout[,4])
a2 <-as.numeric(sortout[,5])
#m1_inf<-as.numeric(sortout[,6])
b<-3^r/2^s
b1<-1*b/b
b2<-0.5*b/b

# 
# 
# #fig2,600by500
# #plots of Estimated ratio r/s and bounds
plot(s, a, frame = T, pch = 20, col = "blue", xlab=expression(paste("Minimum Steps: ",s)),
     ylab = expression(paste("Ratio  ",r/s)), type='o', ylim=c(0.3,0.65))
#ylab = expression(paste("Estimated cycle number bounds")), type='o', log='y')
lines(s,a1,col="red")
lines(s,a2,col="red")

#fig3,600by500
#plots of Estimated ratio 3^r/2^s and bounds
plot(s, b, frame = T, pch = 20, col = "blue", xlab=expression(paste("Minimum Steps: ",s)), 
     ylab = expression(paste("Ratio  ",3^r/2^s)), type='o', ylim=c(0.3,1))
#ylab = expression(paste("Estimated cycle number bounds")), type='o', log='y')
lines(s,b1,col="red")
lines(s,b2,col="red")


i=which.max(a)
c(s[i],r[i],a[i])

#reduce the numbers for Table 4
sortout<-out1[order(out[,3],decreasing=T),]

library(dplyr)

df <- as.data.frame(sortout)
df$err<-log10(df$a2-df$r/df$s)

col_index <- 1
rows_to_delete<-99999999  #set a fake row to be removed


while (length(rows_to_delete)>=1 || rows_to_delete!=0 ){
  # Find the rows where the value is greater than the previous row
  rows_to_delete <- which(df[[col_index]] > lag(df[[col_index]], default = first(df[[col_index]])))
  if (length(rows_to_delete)==0) {break} 
  # Delete those rows
  df <- df[-rows_to_delete, ]
}


# Table 4
options(digits=22)
Tab4<-df[order(df$s), ]



#6/25/2025
#method 2 r1 (very fast version)  to get numbers (r,s) of sequence, get Table 4 numbers
#	Mike Winkler, On the structure and the behaviour of Collatz 3n + 1 sequences, arXiv:1412.0519v1 [math.GM] 27 Nov 2014
# avoid calculating 3^r/2^s
# using formula: s=floor(1+r/a2)

#r range: 
r1=1
r2=400000000 #about 24 min
r2=40000 #about 24 min

out <- matrix(NA, nrow=r2-r1+1, ncol=3)
colnames(out) <- c('s','r','r/s')

tab4 <- matrix(NA, nrow=50, ncol=3)
colnames(tab4) <- c('s','r','r/s')
j=1
tab4[j,]<-cbind(2,1,0.5)


r=0
s=0
a2=mpfr(log(2)/log(3),prec=200)
a2=log(2)/log(3)

start_time <- Sys.time()
for (i in r1:r2){
  r=i
  s=floor(1+r/a2)
  a=r/s
 # a1=(1-1/s)*a2
  #  m1_inf=(3^(r-1)-2^(r-1))/(2^s-3^r)
  out[i-r1+1,]<-cbind(s,r,a)
  
  #find monotonally increase in r/s
  if (a>tab4[j,3]) {
    j=j+1
    tab4[j,]<-cbind(s,r,a)
  }
}

end_time <- Sys.time()
end_time-start_time

tab41<-cbind(tab4,(1-1/tab4[,1])*a2,a2, log10(a2-tab4[,2]/tab4[,1]))
colnames(tab41) <- c('s','r','r/s','a1','a2','d')

write.table(tab41, file="C:/Users/Owner/OneDrive/3n+1/Tab41.csv", sep = ',', row.names = F, col.names = T)

options(digits=22)


#7/3/2025, method 2 r2 
#high-precision calculations,

library(Rmpfr)
x <- mpfr(306, precBits = 256)/485; print(x)
x <- mpfr(306, precBits = 256)/485; print(x)

precBits <- 256 
one <- mpfr(1, precBits)
two <- mpfr(2, precBits)
three <- mpfr(3, precBits)
a2 <- log(two) / log(three); print(a2)
cta <- log(three) / log(two); print(cta)


#zero_vals <- mpfr(rep(0, 3*(r2-r1+1)), prec = 50)
#out <- matrix(zero_vals, nrow=r2-r1+1, ncol=3)
#colnames(out) <- c('s','r','r/s')

# storage matrix 
m=3 # total number of columns
# n=41 # total number of rows
 n=5 # total number of rows


zeros <- mpfr(rep(0, m*n), precBits)
tab4 <- matrix(zeros, nrow=n, ncol=m)
colnames(tab4) <- c('s','r','r/s')

b <- mpfr(1, precBits)/2;print(b)
 b <- mpfr(10781274, prec = precBits)/17087915;print(b-a2)
# b <- mpfr(64497107, prec = precBits)/102225496;print(b)
# b <- mpfr(311704261, prec = precBits)/494039565;print(b)

 j=1; tab4[j,]<-cbind(2,1,b)
 j=1; tab4[j,]<-cbind(17087915,10781274,b)
# j=2; tab4[j,]<-cbind(102225496,64497107,b)


#r range: 
r1=1
r1=64497100 #about 24 min
r2=400000000 #about 24 min

start_time <- Sys.time()
start_time 
for (r in r1:r2) {
  if (j==n) {break}
    a3=r/a2
    s=floor(1+a3)
    a=r/s;
    # a1=(1-1/s)*a2
    #  m1_inf=(3^(r-1)-2^(r-1))/(2^s-3^r)
    #out[i-r1+1,]<-cbind(s,r,a)
    
    #find monotonally increase in r/s
    if (a>b & j<n) {
      j=j+1
      tab4[j,]<-cbind(s,r,a)
      b<-a
      delta<-a2-a;
      print(j); print(s); print(r); print(delta,digit=3)
      print(Sys.time()-start_time)
    }
}

end_time <- Sys.time()
end_time 
end_time-start_time

print(tab4, digits = 128)



library(Rmpfr)

# Precision (in bits)
precBits <- 1000

# Define r and s as mpfr (high precision float)
r <- mpfr("431166034846567", precBits)
s <- mpfr("683381996816440", precBits)

# Compute 3^r, 2^s
three_pow_r <- 3 ^ r
two_pow_s <- 2 ^ s

# Compute numerator: s * 3^r
numerator <- s * three_pow_r

# Compute denominator: 2^s - 3^r
denominator <- two_pow_s - three_pow_r

# Final result
result <- numerator / denominator

# Print result in scientific notation
format(result, scientific = TRUE)







# Eliahou, Shalom (1993). "The 3x + 1 problem: new lower bounds on nontrivial cycle lengths". Discrete Mathematics. 118 (1): 45–56. doi:10.1016/0012-365X(93)90052-U.
# Table 2 start with s=24 727, r=15 601 to the end

#maximum of r/s, data also in 3n+1.xlsx file
s=c(
  2,
  5,
  8,
  27,
  46,
  65,
  149,
  233,
  317,
  401,
  485,
   1539,
   2593,
   3647,
   4701,
   5755,
   6809,
   7863,
   8917,
   9971,
   11025,
   12079,
   13133,
   14187,
   15241,
   16295,
   17349,
   18403,
   19457,
   20511,
   21565,
   22619,
   23673,
   24727,
   75235,
   125743,
   301994,
   17087915,
   102225496,
   187363077,
   272500658,
   630138897,
   10439860591,
   114208327604,
   217976794617,
   1411629234715,
   2605281674813,
   3798934114911,
   4992586555009,
   6186238995107,
   7379891435205,
   8573543875303,
   18340740190704,
   101470897268921,
   184601054347138,
   267731211425355,
   350861368503572,
   433991525581789,
   517121682660006,
   600251839738223,
   683381996816440
  # 117972833293231014+69009683580368485,
  # 124207383220472977+72656661496678846,
  # 130441933147714940+76303639412989207,
  # 254649316368187917+148960300909668053
)



r=c(
  1,
  3,
  5,
  17,
  29,
  41,
  94,
  147,
  200,
  253,
  306,
   971,
   1636,
   2301,
   2966,
   3631,
   4296,
   4961,
   5626,
   6291,
   6956,
   7621,
   8286,
   8951,
   9616,
   10281,
   10946,
   11611,
   12276,
   12941,
   13606,
   14271,
   14936,
   15601,
   47468,
   79335,
   190537,
   10781274,
   64497107,
   118212940,
   171928773,
   397573379,
   6586818670,
   72057431991,
   137528045312,
   890638885193,
   1643749725074,
   2396860564955,
   3149971404836,
   3903082244717,
   4656193084598,
   5409303924479,
   11571718688839,
   64021008208555,
   116470297728271,
   168919587247987,
   221368876767703,
   273818166287419,
   326267455807135,
   378716745326851,
   431166034846567
  # 117972833293231014,
  # 124207383220472977,
  # 130441933147714940,
  # 254649316368187917
   )

s-r

c=a2-mpfr(r,prec=256)/s; print(r,digits = 3); print(c,digits = 3)

# Table 2 c=mpfr(s,prec=256)/r - cta; print(r,digits = 3); print(c,digits = 3)


plot(log(r),log(s))
lines(log(r),log(s))

log10(-10781274/17087915+log(2)/log(3))

-171928773/272500658+log(2)/log(3)

sup=s/s*log(2)/log(3)
inf=(1-1/s)*sup

plot(s, c, frame = T, pch = 20, col = "blue", xlab=expression(paste("Minimum Steps: ",s)), 
     ylab = expression(paste("Closest ratio  ",r/s, "  to  ", log(2)/log(3))), type='o',
     xlim=c(0,30500000),
     #ylim=c(0.6309,0.63093),
     #ylim=c(0.08,20000),
     #log='x'
)
lines(s,inf,col="red")
lines(s,sup,col="red")


#fig8,600by500
library(svglite)
svglite("fig8.svg", width = 6, height = 5)

d=log10(c);  print(c,digits = 3)
x<-log10(s)
# d=mpfr(log10(sup-c),prec=256)
# x=mpfr(log10(s),prec=256)


plot(x, d, frame = T, pch = 20, col = "blue", xlab=expression(paste("Minimum Steps (s) in ", 10^x)), 
     #     ylab = expression(paste("Closest distance of  ",r/s, "  to  ", log2/log3)), type='o',
     ylab = expression(paste("Closest distance between r/s and ", log[3](2), " in ", 10^x)), type='o',
     # xlim=c(10,30500000),
     #ylim=c(0.6309,0.63093),
     #ylim=c(1e-13,1e-04),
     # log='y'
)

y <- asNumeric(d)

model <- lm(y ~ x)
# plot(model1)
abline(model, col = "red", lwd = 2)

coefs <- coef(model)
eq <- sprintf("log10(log3(2) - r/s) = %.2f %.2f log10(s)", coefs[1], coefs[2])
r2 <- summary(model)$r.squared
r2_txt <- sprintf("R² = %.3f", r2)

# Add text to plot
text(x = 5, y = max(y), labels = eq, pos = 4, col = "blue")
text(x = 5, y = max(y) - 2, labels = r2_txt, pos = 4, col = "blue")

dev.off() 

new_data <- data.frame(x = c(40, 250, 1000))
predictions <- predict(model, newdata = new_data)

cc <- asNumeric(c)
 fig8dat<-data.frame(s=s,r=r,inf=inf,c=cc,sup=sup, diff=y)

print(cc,digit=3) 


2.7182818284587582*2^3.5*30^5*log(3) 


821013301*log10(683381996816440)
821013301*log2(683381996816440)


# this is computing 1-3^r/2^s
1-exp(r*log(three)-s*log(two))

r <- mpfr("148960300909668053", precBits)
s <- mpfr("254649316368187917", precBits)
r/s

r <- mpfr("431166034846567", precBits)
s <- mpfr("683381996816440", precBits)
1-exp(r*log(three) - s*log(two) )

 
# min m1:
r=306
s=485

min_m1 >= ( 3^(r) - 2^(r) )/(2^s - 3^r) 
min_m1 > ( 3^(r) - 1/(2^(1-s*log(2)/log(3) ) +s)) / (2^s - 3^r) 
# 977.7447725446728554743

# here , we get :
# min_m1>=326

eN=(1-3^r/2^s)

min_m1 = 1/eN

1 > (3^(r) - 2^(r))/2^s



# Hercher, C. (2023). "There are no Collatz m-cycles with m <= 91" (PDF). Journal of Integer Sequences. 26 (3): Article 23.3.5.

# Let m be the number of local minima in such a nontrivial cycle. Simons and de
# Weger proved that m ≥ 76. With newer bounds on the range of starting values for
# which the Collatz conjecture has been checked, one gets m ≥ 83. In this paper, we
# prove m ≥ 92.





plot(s, a, frame = T, pch = 20, col = "blue", xlab=expression(paste("Minimum Steps: ",s)), 
     ylab = expression(paste("Closest ratio  ",r/s, "  to  ", log(2)/log(3))), type='o',
            xlim=c(0,30500000),
     #ylim=c(0.6309,0.63093),
     #ylim=c(0.08,20000),
     #log='x'
    )
lines(s,inf,col="red")
lines(s,sup,col="red")

olddata<- data.frame(s=s, r=r)

#method 3: using linear reg results to get numbers (r,s) of sequence, get Table 4 numbers
# avoid 3^r/2^s

# linear regreesion of s on r, 3/30/2024
# may generate more (s,r) pairs

model1 <- lm(s ~ r)
# plot(model1)

round(predict(model1,newdata=list(s=NA,r=190537)))

newdata<-data.frame(s=c(s,NA),r=c(r,190537))

newdata$r/round(predict(model1, newdata = newdata))

log(2)/log(3)


summary(model1)

r=190537
# method 1:
  s_hat <- round(6.690405800466725636e-04+1.584962500685801690e+00*r)
# method 2:
  s_hat <- round(predict(model1,newdata=list(s=NA,r=r)))

  start_time <- Sys.time()
  rmax <-2000000
  out <- matrix(NA, nrow=rmax, ncol=5)
  colnames(out) <- c('s','r','r/s','a1','a2')

  a2=log(2)/log(3)  
  r0=100*190537+1

  s0=floor(r/a2)
r/s

  for (i in r0:(r0+rmax) ){
    r=i
    # s=r
    s=round(9.201365614150338202e-04+1.584962482180414245e+00*r)-1
    a=r/s
    a1=(1-1/s)*a2
    
    while (a<a1 | a>a2)
    {
      a=r/s
      a1=(1-1/s)*a2
      s=s+1
    }
    s<-s-1
    a=r/s
    a1=(1-1/s)*a2
    #  m1_inf=(3^(r-1)-2^(r-1))/(2^s-3^r)
    out[(i-r0),]<-cbind(s,r,a,a1,a2)
  }
  
  #c=0.630927835
  #c=0.63092975
  c=0

  out1<-out[which(out[,3]>=c),]
  sortout<-out1[order(out1[,1],decreasing=F),]
  
  s <-as.numeric(sortout[,1])
  r <-as.numeric(sortout[,2])
  a <-as.numeric(sortout[,3])

  i=which.max(a)
  c(s[i],r[i],a[i])

  190537/301994
  # [1] 0.6309297535712630189408
  19164902/30375651 
  # [1] 0.6309297535713720428419
      # 0.6309297535714574189925
  end_time <- Sys.time()
  end_time - start_time
  
# pred<- exp(-0.4115124337836175705+1.6350958356783278358*p2)
# lines(p2,pred,col="green") #Add a regression line

conf_interval <- exp(predict(model1, newdata=out, interval="confidence",level = 0.95))
lines(Fm, conf_interval[,1], col="red", type = "l",lty=1,lwd=4)






exp(-5*(log(2)^5))



#fig1:600by500

library(svglite)
svglite("fig1.svg", width = 6, height = 5)

# With ggsave()
#ggsave("myplot.svg", width = 8, height = 8, units = "cm")


s=c(1,2,3,4,5,6,7,8)
s=c(0,1,2,3,4,5,6,7)
m=c(187,	281,	422,	211,	317,	476,	238,	119)
c1=m[1]*m/m 
c2=c1/2
plot(s, m, frame = T, pch = 20, col = "blue", xlab=expression(paste("Step (s)")), 
      ylab = expression(paste(" ",Fq(m))), type='o', ylim=c(0,500))
 abline(h=m[1], col="red")
 abline(h=m[1]/2, col="red")
 text(6.6, 238, '2Fq(m)')
 text(6.5, 119, 'Fq(m)')
 text(0.5, 200, 'm') 
 text(0.5, 110, 'm/2')  
 text(5, 205, 'q=1101100')  
 
dev.off() 
 
N=485

eN=2^(-(1-log(2)/log(3))*N-1)
e=1-0.99999999999999897828317652204
M=1/3/e-1/3-eN/e


C<-2.7182818284587582*2^3.5*30^5*log(3) 
B=1/(2.7182818284587582*485)^1.06

log(2)



# in ref [6]
B=s/log(3)+r
H=max(log(B)+0.18,10)

H=log(B)+0.18
d=2^(-32.31*H^2*log2(3)-1)



1/2^(s-r+1)


# 4/18/2024

# Eliahou, Shalom (1993). "The 3x + 1 problem: new lower bounds on nontrivial cycle lengths". Discrete Mathematics. 118 (1): 45–56. doi:10.1016/0012-365X(93)90052-U.
# Table1 from equation: (3.2)

a<-c(1,1,1,2,2,3,1,5,2,23,2,2,1,1,55,1,4,3,1,1,15,1,9,2,5,7,1,1,4,8,1,11,1,20,2,1)
length(a)
p<-0
q<-0
p[1]<-0
p[2]<-1
q[1]<-1
q[2]<-0

for (n in 3:38){
  p[n]<-a[n-2]*p[n-1]+p[n-2]
  q[n]<-a[n-2]*q[n-1]+q[n-2]
}

p/q-log(3)/log(2)

# max numbers
log10(log(2)/log(3)-10781274/17087915)
log10(log(2)/log(3)-171928773/272500658)

(1-1/17087915)*log(2)/log(3)
(1-1/272500658)*log(2)/log(3)



111202/176251
397513379/630138897

53715833/85137581-log(2)/log(3)
6586818670/10439860591-log(2)/log(3)
137528045312/217976794617-log(2)/log(3)

267118416222671843/423372672964960618-log(2)/log(3)
397560349370386783/630118245525664765-log(2)/log(3)


# Table 2

n<-c(9,9,11,13,15,15,15,17,19,21,21,23,23,23,23,23,23,23,25,27,27,27,27,27,27,27,27,29)
i<-c(0,1,0, 0, 0, 1, 2, 0, 0, 0, 1, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 2, 3, 4, 5, 6, 7, 0)
length(n)
length(i)


p_ni<-0
q_ni<-0
# p_ni[1]<-24727
# q_ni[1]<-15601
# p_ni[2]<-125743
# q_ni[2]<-79335

for (j in 1:28){
  k=n[j]+3
  l=i[j]
  p_ni[j]<-p[k]+l*p[k+1]
  q_ni[j]<-q[k]+l*q[k+1]
}

p_ni/q_ni-log(3)/log(2)


r<-c(2,
     3,
     4,
     5,
     6,
     7,
     8,
     9,
     10,
     11,
     12,
     13,
     14,
     15,
     16,
     17,
     18,
     19,
     20,
     21,
     22,
     23,
     24,
     25,
     26,
     27,
     28,
     29,
     30,
     306)
s<-c(4,
     5,
     7,
     8,
     10,
     12,
     13,
     15,
     16,
     18,
     20,
     21,
     23,
     24,
     26,
     27,
     29,
     31,
     32,
     34,
     35,
     37,
     39,
     40,
     42,
     43, 
     45,
     46,
     48,
     485)


min_m1=10*(3^r-2^r)/(2^s-3^r)
max_m1 = s*3^r/(2^s-3^r)  # Ref [7] Florian Luca, On the non-trivial cycles in Collatz’s problem, SUT Journal of Mathematics Vol. 41, No. 1 (2005), 31–41

max_m1=(3^(r-2)*5+(r-2)*2^(s-2) )/(2^s-3^r)


C<-2.7182818284587582*2^3.5*30^5*log(3) 
B=1/(2.7182818284587582*485)^1.06
e=2.7182818284587582
e*17087915
log2(log2(e*17087915)*C) 
log10(log10(e*17087915)*C)

