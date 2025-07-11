# -----------------------------------------------------------
# Program: theoretical accuracy
# Details:
#   Functions to calculate the theoretical accuracy 
#   for each combination of decision strategy and evaluation criterion
#  
#   DDEC_S: Deterministic decision (DD)’s accuracy under 
#           evaluation criterion based on the sub-population mean (EC-S)
#   RDEC_S: Randomized decision (RD)’s accuracy under EC-S 
#   DDEC_P: DD’s accuracy under evaluation criterion based on the population mean (EC-P)
#   RDEC_P: RD’s accuracy under EC-P
#
#
# Input Arguments:
#   mu  : population mean difference
#         (u1 - u2; mean of population 1 minus mean of population 2)
#
#   s1  : standard deviation of X1
#   s2  : standard deviation of X2 
#         X_t=g(X_t1,...,X_(tn_t)), where t=1,2. 
#         Only two types of function g(.) are considered: 
#         (1) selecting a particular random sample in the set {X_t1,...,X_(tn_t)}
#         s1=standard deviation of sub-population 1 
#         s2=standard deviation of sub-population 2 
#         (2) selecting the mean, X_t=Xbar_t.   
#         s1=standard deviation of sub-population 1 divided by square root of n1 
#         s2=standard deviation of sub-population 2 divided by square root of n2
#        
#   r1  : standard deviation of population 1
#   r2  : standard deviation of population 2
#
# Output:
#   Returns a list containing:
#     AC : accuracy 
# -----------------------------------------------------------

# Load libraries 
library(cubature)
library(tidyr)

# -----------------------------------------------------------
# Deterministic decision (DD)’s accuracy under 
# evaluation criterion based on the sub-population mean (EC-S)
# -----------------------------------------------------------
DDEC_S<-function(mu,s1,s2,r1,r2){
  # AC_ds=P(X>0,m>0)+P(X<=0,m<=0)

  #f_X,m(x, w)=f_X|m(x|w)f_m(w)
  # x[1]=X=X1-X2
  # x[2]=m=m1-m2
  f_xm = function(x){
    (1/(2*pi*sqrt(s1^2+s2^2)*sqrt(r1^2+r2^2))) * exp( -0.5* ( (x[1] - x[2])/sqrt(s1^2+s2^2) )^2 - 0.5* ( ((x[2]-mu)^2)/(r1^2+r2^2) ))
  }
  
  #P(X>0,m>0)
  case1 = adaptIntegrate(f_xm, lowerLimit = c(0, 0), upperLimit = c(Inf, Inf))
  #P(X<=0,m<=0)
  case2 = adaptIntegrate(f_xm, lowerLimit = c(-Inf, -Inf), upperLimit = c(0, 0))
  
  # AC_ds
  rs<-case1$integral + case2$integral 
  return(list(AC=rs))
}

# -----------------------------------------------------------
# Randomized decision (RD)’s accuracy under EC-S
# -----------------------------------------------------------
RDEC_S<-function(mu,s1,s2,r1,r2){ 
  #AC_rs=first term + second term in paper
  #first term
  RDECS1<-function(a){
    # f_m|X(w|x)
    # a=x=x1-x2
    # x[1]=m=m1-m2
    f1 = function(x){
      (1/sqrt(2*pi*(((s1^2+s2^2)*(r1^2+r2^2))/(r1^2+s1^2+r2^2+s2^2)))) * exp( -0.5* (1/(((s1^2+s2^2)*(r1^2+r2^2))/(r1^2+s1^2+r2^2+s2^2)))*( (x[1]- (((r1^2+r2^2)/(r1^2+s1^2+r2^2+s2^2))*a+((s1^2+s2^2)/(r1^2+s1^2+r2^2+s2^2))*mu))^2) ) 
    }
    #f_X,m(x,w)=f_X|m(x|w)f_m(w)
    # a=x=x1-x2
    # x[1]=m=m1-m2
    f2 = function(x){
      (1/(2*pi*sqrt(s1^2+s2^2)*sqrt(r1^2+r2^2))) * exp( -0.5* ( (a - x[1])/sqrt(s1^2+s2^2) )^2 - 0.5* ( ((x[1]-mu)^2)/(r1^2+r2^2) )  )
    }
    
    # integrate from w=0 to Inf
    case1 = adaptIntegrate(f1, lowerLimit = c(0), upperLimit = c(Inf) )
    
    case11 = adaptIntegrate(f2, lowerLimit = c(0), upperLimit = c(Inf) )
    
    rs<-case1$integral *case11$integral 
    return(rs)
  }
  
  #second term
  RDECS2<-function(a){
    f1 = function(x){
      (1/sqrt(2*pi*(((s1^2+s2^2)*(r1^2+r2^2))/(r1^2+s1^2+r2^2+s2^2)))) * exp( -0.5* (1/(((s1^2+s2^2)*(r1^2+r2^2))/(r1^2+s1^2+r2^2+s2^2)))*( (x[1]-(((r1^2+r2^2)/(r1^2+s1^2+r2^2+s2^2))*a+((s1^2+s2^2)/(r1^2+s1^2+r2^2+s2^2))*mu))^2))
    }
    
    f2 = function(x){
      (1/(2*pi*sqrt(s1^2+s2^2)*sqrt(r1^2+r2^2))) * exp( -0.5* ( (a - x[1])/sqrt(s1^2+s2^2) )^2 - 0.5* ( ((x[1]-mu)^2)/(r1^2+r2^2) )  )
    }
    
    # integrate from w=-Inf to 0 
    case1 = adaptIntegrate(f1, lowerLimit = c(-Inf), upperLimit = c(0) )
    
    case11 = adaptIntegrate(f2, lowerLimit = c(-Inf), upperLimit = c(0) )
    
    rs<-case1$integral *case11$integral 
    return(rs)
  }
  
  # integrate from x=-Inf to Inf  
  case1pre = adaptIntegrate(RDECS1, lowerLimit = c(-Inf), upperLimit = c(Inf) )
  
  case2pre = adaptIntegrate(RDECS2, lowerLimit = c(-Inf), upperLimit = c(Inf) )
  
  case1<-case1pre$integral
  case2<-case2pre$integral
  
  # AC_rs
  rs<-case1+case2
  return(list(AC=rs))
} 

# -----------------------------------------------------------
# DD’s accuracy under 
# Evaluation criterion based on the population mean (EC-P)
# -----------------------------------------------------------

DDEC_P<-function(mu,s1,s2,r1,r2){
  # AC_dp=P(X>0,mu>0)+P(X<=0,mu<=0)
  DDECP1<-function(a){
  #f_X,m(x, w)=f_X|m(x|w)f_m(w)
  # mu=u=u1-u2
  # a=X=X1-X2
  # x[1]=m=m1-m2
  f3 = function(x){
    (1/(2*pi*sqrt(s1^2+s2^2)*sqrt(r1^2+r2^2))) * exp( -0.5* ( (a - x[1])/sqrt(s1^2+s2^2) )^2 - 0.5* ( ((x[1]-mu)^2)/(r1^2+r2^2)))
  }
  
  # integrate from w=-Inf to Inf 
  case1 = adaptIntegrate(f3, lowerLimit = c(-Inf), upperLimit = c(Inf) )
  
  rs<-case1$integral 
  return(rs)
  }

  if (mu>0){ 
    #P(X>0)
    case = adaptIntegrate(DDECP1, lowerLimit = c(0), upperLimit = c(Inf) )
  } else{ 
    #P(X<=0)
    case = adaptIntegrate(DDECP1, lowerLimit = c(-Inf), upperLimit = c(0) )
  }
  
  # AC_dp
  rs<-case$integral 
  return(list(AC=rs))
}

# -----------------------------------------------------------
# RD’s accuracy under EC-P
# -----------------------------------------------------------

RDEC_P<-function(mu,s1,s2,r1,r2){ 

  RDECP1<-function(a){  
    # f_m|X(w|x)
    # a=x=x1-x2
    # x[1]=m=m1-m2
    f1 = function(x){
      (1/sqrt(2*pi*(((s1^2+s2^2)*(r1^2+r2^2))/(r1^2+s1^2+r2^2+s2^2)))) * exp( -0.5* (1/(((s1^2+s2^2)*(r1^2+r2^2))/(r1^2+s1^2+r2^2+s2^2)))*( (x[1]- (((r1^2+r2^2)/(r1^2+s1^2+r2^2+s2^2))*a+((s1^2+s2^2)/(r1^2+s1^2+r2^2+s2^2))*mu))^2) ) 
    }
    #f_X,m(x,w)=f_X|m(x|w)f_m(w)
    # a=x=x1-x2
    # x[1]=m=m1-m2
    f111 = function(x){
      (1/(2*pi*sqrt(s1^2+s2^2)*sqrt(r1^2+r2^2))) * exp( -0.5* ( (a - x[1])/sqrt(s1^2+s2^2) )^2 - 0.5* ( ((x[1]-mu)^2)/(r1^2+r2^2) )  )
    }
    #f_X(x)
    fx= adaptIntegrate(f111, lowerLimit = c(-Inf), upperLimit = c(Inf) )
    
    if (mu>0){
      #P(m>0|x)
      case1 = adaptIntegrate(f1, lowerLimit = c(0), upperLimit = c(Inf) )
      
    }
    if (mu<=0){
      #P(m<=0|x)
      case1 = adaptIntegrate(f1, lowerLimit = c(-Inf), upperLimit = c(0) )
    }
    rs1<-case1$integral*fx$integral
    return(rs1)
  }
  
  # integrate from x=-Inf to Inf
  casepre = adaptIntegrate(RDECP1, lowerLimit = c(-Inf), upperLimit = c(Inf) )
  
  #AC_rp
  rs<-casepre$integral
  return(list(AC=rs))
}


# -----------------------------------------------------------
# Example 
# -----------------------------------------------------------

dd.s = DDEC_S(0.5,1,1,1,1)
rd.s = RDEC_S(0.5,1,1,1,1)
dd.p = DDEC_P(0.5,1,1,1,1)
rd.p = RDEC_P(0.5,1,1,1,1)

dd.s
rd.s
dd.p
rd.p
