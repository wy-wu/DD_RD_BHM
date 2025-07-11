# -----------------------------------------------------------
# Program: acc.calculate.r
# Details:
#         Section 4 – Simulated Case Studies:
#         4.1 Computes accuracy under the assumption that the true parameter values are known.
#         4.2 Performs decision-making based on observed data and estimated parameters.
#             Note: Three scenarios are evaluated; please refer to the manuscript for details.
# -----------------------------------------------------------


library(cubature)
source("theory_AC.R")

# -----------------------------------------------------------
# 4.1 Calculating the accuracy assuming that the true parameter
#     values are known
# -----------------------------------------------------------

dd.s = DDEC_S(-0.04,0.07,0.08,0.05,0.10)
rd.s = RDEC_S(-0.04,0.07,0.08,0.05,0.10)
dd.p = DDEC_P(-0.04,0.07,0.08,0.05,0.10)
rd.p = RDEC_P(-0.04,0.07,0.08,0.05,0.10)
dd.s
rd.s
dd.p
rd.p



# -----------------------------------------------------------
# 4.2 Decision making based on the observed data and 
#     estimated parameter values
# -----------------------------------------------------------



### Define function for eq5
u1=0.23; u2=0.27; mu = u1-u2
r1=0.05; r2=0.10
s1=0.07; s2=0.08 
x1=0.24; x2=0.26; a=x1-x2

f1 = function(x){
  (1/sqrt(2*pi*(((s1^2+s2^2)*(r1^2+r2^2))/(r1^2+s1^2+r2^2+s2^2)))) * 
   exp( -0.5* (1/(((s1^2+s2^2)*(r1^2+r2^2))/(r1^2+s1^2+r2^2+s2^2)))*
 ( (x[1]- (((r1^2+r2^2)/(r1^2+s1^2+r2^2+s2^2))*a+((s1^2+s2^2)/(r1^2+s1^2+r2^2+s2^2))*mu))^2) ) 
}

# integrate from w=0 to Inf
case1 = adaptIntegrate(f1, lowerLimit = c(0), upperLimit = c(Inf) )
case1

#############
# Scenario 1
#############
x1=0.24; x2=0.26; a=x1-x2
case1 = adaptIntegrate(f1, lowerLimit = c(0), upperLimit = c(Inf) )
case1

#############
# Scenario 2
#############
x1=0.25; x2=0.24; a=x1-x2
case2 = adaptIntegrate(f1, lowerLimit = c(0), upperLimit = c(Inf) )
case2


##############
# Scenario 3
##############
## From example3.R, we have the estimated parameters:
# u_1=0.24,u2=0.26, r1=0.03, r2=0.08,s1=0.03 and s2=0.08.

x1 = mean(c(0.23, 0.25, 0.28))
x1
mean(c(0.19, 0.22, 0.23, 0.27))
mean(c(0.21, 0.22, 0.24, 0.28))

mean(c(0.20, 0.29))
x2 = mean(c(0.19, 0.28, 0.34))  
x2

# since x1 and x2 are mean
# s1= s1_hat/sqrt(3)
# s2= s2_hat/sqrt(3)
u1=0.24; u2=0.26; mu = u1-u2
r1=0.03; r2=0.08
s1=0.03/sqrt(3); s2=0.08/sqrt(3) 
x1=0.253; x2=0.27; a=x1-x2

case3a = adaptIntegrate(f1, lowerLimit = c(0), upperLimit = c(Inf) )
case3a

###
# s1= s1_hat/sqrt(3)
# s2= s2_hat/sqrt(2)

u1=0.24; u2=0.26; mu = u1-u2
r1=0.03; r2=0.08
s1=0.03/sqrt(3); s2=0.08/sqrt(2) 
x1=0.253; x2=0.245; a=x1-x2

case3b = adaptIntegrate(f1, lowerLimit = c(0), upperLimit = c(Inf) )
case3b


