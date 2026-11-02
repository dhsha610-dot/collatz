# environment and functions for Collatz process 
options(digits=22)

library(gmp)  #for big numbers
library(Rmpfr) #for small numbers
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


is.even <- function(x) {
  if (as.bigz(x) %% 2 == 0) TRUE
  else FALSE
}

collatz <- function(n) {
  if (is.even(n)) as.bigz(n)/2
  else (3 * as.bigz(n) + 1)/2
}

bigz1=as.bigz("27353174180740493936925409695987002025536604159003070399660429096649393310061075405506794475763653718625318901355796621682534492976519463273944176643274761783713821522632289899076598818508786143223704385779310297725363258434632088907939839")
bigz2=as.bigz("37156274889223673328946970673189379090504259790996737861419817121864723982393460680622326589389941666295331379135808502905715522514535613969322698413735982793541065938680832420603238792168680545393321856628980319639419302924341360753269755") 

# get number of digits for a big z
nchar(as.character(bigz1))
nchar(as.character(bigz2))

is.even(bigz1)
collatz(bigz1)

is.even(bigz2)
collatz(bigz2)

big_1=as.bigz(1)
big_2=as.bigz(2)

big_number = as.bigz("3478481451834367")

