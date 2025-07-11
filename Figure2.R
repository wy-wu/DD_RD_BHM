# -----------------------------------------------------------
# Program: Figure2
# Details:
#   The theoretical (curves) and simulation (dots) accuracy:
#   AC_ds: Deterministic decision (DD)'s accuracy under evaluation criterion based on the sub-population mean (EC-S)
#   AC_rs: Randomized decision (RD)'s accuracy under evaluation criterion based on the sub-population mean (EC-S) 
#   AC_dp: DD's accuracy under EC-S
#   AC_rp: RD's accuracy under EC-P 
#
#   We explore five distinct scenarios characterized by varying 
#   population variances(r_1 and r_2) and sub-population variances (s_1 and s_2).
#   
#   Scenario (variance): 
#   1. Population = Sub-population;	                r_1=r_2=1	s_1=s_2=1
#   2. Population same, Sub-population different; 	r_1=r_2=1	s_1=1,s_2=2
#   3. Population < Sub-population;	                r_1=r_2=1	s_1=s_2=2
#   4. Population different, Sub-population same;	  r_1=1,r_2=2	s_1=s_2=1
#   5. Population > Sub-population; 	              r_1=r_2=2	s_1=s_2=1
#    
# -----------------------------------------------------------

# Load libraries 
library(cubature)
library(tidyr)
library(ggplot2)
library(patchwork)

source("theory_AC.R")

# -----------------------------------------------------------
# Scenario1 
# -----------------------------------------------------------

# compute the theoretical and simulation accuracy

compute_accuracy_data <- function(
    u1 = 7,
    s1 = 1,
    s2 = 1,
    r1 = 1,
    r2 = 1
) {
  # Theoretical accuracy
  # Define the grid of u2 
  v_u2_thry <-sort(c(seq(1,13,by=1),seq(0, 14, length=1000)))
  # Loop through the values:
  oc_list_thry <- lapply(v_u2_thry, function(x) {
    dd.s <- DDEC_S(u1 - x, s1, s2, r1, r2)
    rd.s <- RDEC_S(u1 - x, s1, s2, r1, r2)
    dd.p <- DDEC_P(u1 - x, s1, s2, r1, r2)
    rd.p <- RDEC_P(u1 - x, s1, s2, r1, r2)
    data.frame(
      u2 = x,
      AC_ds = dd.s$AC,
      AC_rs = rd.s$AC,
      AC_dp = dd.p$AC,
      AC_rp = rd.p$AC
    )
  })
  oc_thry <- do.call(rbind, oc_list_thry)
  
  # Reshape theoretical data
  oc_thry_EC.S <- tidyr::gather(oc_thry, "Strategy and EC ", "AC", c(2,3))
  oc_thry_EC.P <- tidyr::gather(oc_thry, "Strategy and EC ", "AC", c(4,5))
  
  # Simulation accuracy
  # ntr: number of iteration 
  simu<-function(u1,u2,s1,s2,r1,r2,ntr,seed){
    # simulate X1,X2,m1,m2  
    if(!is.na(seed)) set.seed(seed)
    m1<-rep(NA,ntr)
    m2<-rep(NA,ntr)
    X1<-rep(NA,ntr)
    X2<-rep(NA,ntr)
    Prob_t<-rep(NA,ntr)
    
    for(i in 1:ntr){
      # generate one value of m1, given u1 and r1.
      # generate one value of m2, given u2 and r2.
      m1[i]<-rnorm(1,u1,r1)
      m2[i]<-rnorm(1,u2,r2)
      #Given m1 and s1,generate  one x1.
      #Given m2 and s2, generate one x2.
      X1[i]<-rnorm(1,m1[i],s1)
      X2[i]<-rnorm(1,m2[i],s2)
      # Compute prob(m1>m2)=prob(m1-m2>0) based on the posterior distribution 
      #Theoretically
      Prob_t[i]<- 1-pnorm(0, ((r1^2+r2^2)/(s1^2+r1^2+s2^2+r2^2))*(X1[i]-X2[i])+((s1^2+s2^2)/(s1^2+r1^2+s2^2+r2^2))*(u1-u2), sqrt(((r1^2+r2^2)*(s1^2+s2^2))/(s1^2+r1^2+s2^2+r2^2)))
    }
    
    output<-cbind(X1,X2,m1,m2,Prob_t)
    output<-as.data.frame(output)
    
    #DD
    output$DD<-ifelse(output$X1>output$X2,"t=1","t=2")
    #RD 
    output$RD<-ifelse(sapply(output$Prob_t, function(x) rbinom(1, 1, x))==1,"t=1","t=2")
    
    
    #EC-S
    output$EC.s<-ifelse(output$m1>output$m2,"t=1","t=2")
    #EC-P
    output$EC.p<-ifelse(u1>u2,"t=1","t=2")
    
    #AC_ds
    AC_ds<-mean(output$DD==output$EC.s)
    #AC_dp
    AC_dp<-mean(output$DD==output$EC.p)
    #AC_rs
    AC_rs<-mean(output$RD==output$EC.s)
    #AC_rp
    AC_rp<-mean(output$RD==output$EC.p)
    
    rs<-data.frame(u2=u2,AC_ds=AC_ds,AC_rs=AC_rs,AC_dp=AC_dp,AC_rp=AC_rp)
    return(rs)
  }
  
  
  v_u2_simu <-seq(0, 14, by=1)
  oc_list_simu <- lapply(v_u2_simu, function(x) {
    simu(
      u1 = u1,
      u2 = x,
      s1 = s1,
      s2 = s2,
      r1 = r1,
      r2 = r2,
      ntr = 100000,
      seed = 1
    )
  })
  oc_simu <- do.call(rbind, oc_list_simu)
  
  # Reshape simulation data
  oc_simu_EC.S <- tidyr::gather(oc_simu, "Strategy and EC ", "AC", c(2,3))
  oc_simu_EC.P <- tidyr::gather(oc_simu, "Strategy and EC ", "AC", c(4,5))
  
  # Return as a list
  return(list(
    oc_thry_EC.S = oc_thry_EC.S[,-c(2,3)],
    oc_thry_EC.P = oc_thry_EC.P[,-c(2,3)],
    oc_simu_EC.S = oc_simu_EC.S[,-c(2,3)],
    oc_simu_EC.P = oc_simu_EC.P[,-c(2,3)]
  ))
}


S1_data<-compute_accuracy_data(u1 = 7,
                               s1 = 1,
                               s2 = 1,
                               r1 = 1,
                               r2 = 1)


# ready to plot
plot_accuracy <- function(
    data_thry,
    data_simu,
    colors = c("#FF00FF", "cyan3"),
    fill_colors = c("#FF00FF", "cyan3"),
    y_lab="Accuracy"
){ggplot() + 
    geom_vline(xintercept=7, linetype=1,linewidth=0.5, colour="lightgray") +
    geom_line(data = data_thry, size=0.8,aes(x=u2, y=AC,colour=`Strategy and EC `, linetype=`Strategy and EC `)) +
    geom_point(data = data_simu, 
               aes(x = u2, y = AC, fill=`Strategy and EC `,shape = `Strategy and EC `))+
    
    scale_shape_manual(values=c("circle filled","circle filled"),name="")+
    scale_linetype_manual(values=c("solid","solid"), name="Strategy and EC")+
    scale_colour_manual(values=colors, name="Strategy and EC") + 
    scale_fill_manual(values = fill_colors)+
    guides(shape = "none", fill = "none")  + 
    scale_y_continuous(breaks=seq(0, 1, by=0.05),limits=c(0.5,1)) + 
    scale_x_continuous(breaks=seq(0, 14, by=2))+
    labs(x=expression(mu[2]), y=y_lab, title="") + 
    theme_bw(base_size=12)}

#AC_ds,AC_rs 
s1_EC.S <- plot_accuracy(
  data_thry = S1_data$oc_thry_EC.S,
  data_simu = S1_data$oc_simu_EC.S,
  colors = c("#FF00FF", "cyan3"),
  fill_colors = c("#FF00FF", "cyan3"),
  y_lab="Accuracy"
)

# AC_dp,AC_rp 
s1_EC.p <- plot_accuracy(
  data_thry = S1_data$oc_thry_EC.P,
  data_simu = S1_data$oc_simu_EC.P,
  colors = c("red", "navyblue"),
  fill_colors = c("red", "navyblue"),
  y_lab="Accuracy"
)


# -----------------------------------------------------------
# Scenario2 
# -----------------------------------------------------------

S2_data<-compute_accuracy_data(u1 = 7,
                               s1 = 1,
                               s2 = 2,
                               r1 = 1,
                               r2 = 1)


# ready to plot
#AC_ds,AC_rs 
s2_EC.S <- plot_accuracy(
  data_thry = S2_data$oc_thry_EC.S,
  data_simu = S2_data$oc_simu_EC.S,
  colors = c("#FF00FF", "cyan3"),
  fill_colors = c("#FF00FF", "cyan3"),
  y_lab=NULL
)

# AC_dp,AC_rp 
s2_EC.p <- plot_accuracy(
  data_thry = S2_data$oc_thry_EC.P,
  data_simu = S2_data$oc_simu_EC.P,
  colors = c("red", "navyblue"),
  fill_colors = c("red", "navyblue"),
  y_lab=NULL
)
# -----------------------------------------------------------
# Scenario3 
# -----------------------------------------------------------
S3_data<-compute_accuracy_data(u1 = 7,
                               s1 = 2,
                               s2 = 2,
                               r1 = 1,
                               r2 = 1)


# ready to plot
#AC_ds,AC_rs 
s3_EC.S <- plot_accuracy(
  data_thry = S3_data$oc_thry_EC.S,
  data_simu = S3_data$oc_simu_EC.S,
  colors = c("#FF00FF", "cyan3"),
  fill_colors = c("#FF00FF", "cyan3"),
  y_lab=NULL
)

# AC_dp,AC_rp 
s3_EC.p <- plot_accuracy(
  data_thry = S3_data$oc_thry_EC.P,
  data_simu = S3_data$oc_simu_EC.P,
  colors = c("red", "navyblue"),
  fill_colors = c("red", "navyblue"),
  y_lab=NULL
)

# -----------------------------------------------------------
# Scenario4 
# -----------------------------------------------------------
S4_data<-compute_accuracy_data(u1 = 7,
                               s1 = 1,
                               s2 = 1,
                               r1 = 1,
                               r2 = 2)


# ready to plot
#AC_ds,AC_rs 
s4_EC.S <- plot_accuracy(
  data_thry = S4_data$oc_thry_EC.S,
  data_simu = S4_data$oc_simu_EC.S,
  colors = c("#FF00FF", "cyan3"),
  fill_colors = c("#FF00FF", "cyan3"),
  y_lab=NULL
)

# AC_dp,AC_rp 
s4_EC.p <- plot_accuracy(
  data_thry = S4_data$oc_thry_EC.P,
  data_simu = S4_data$oc_simu_EC.P,
  colors = c("red", "navyblue"),
  fill_colors = c("red", "navyblue"),
  y_lab=NULL
)


# -----------------------------------------------------------
# Scenario5 
# -----------------------------------------------------------
S5_data<-compute_accuracy_data(u1 = 7,
                               s1 = 1,
                               s2 = 1,
                               r1 = 2,
                               r2 = 2)


# ready to plot
#AC_ds,AC_rs 
s5_EC.S <- plot_accuracy(
  data_thry = S5_data$oc_thry_EC.S,
  data_simu = S5_data$oc_simu_EC.S,
  colors = c("#FF00FF", "cyan3"),
  fill_colors = c("#FF00FF", "cyan3"),
  y_lab=NULL
)

# AC_dp,AC_rp 
s5_EC.p <- plot_accuracy(
  data_thry = S5_data$oc_thry_EC.P,
  data_simu = S5_data$oc_simu_EC.P,
  colors = c("red", "navyblue"),
  fill_colors = c("red", "navyblue"),
  y_lab=NULL
)



# -----------------------------------------------------------
# Combine all plots into a single figure
#   Top row: EC-S plots for Scenarios 1–5
#   Bottom row: EC-P plots for Scenarios 1–5
# -----------------------------------------------------------


# Define common theme 
my_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size=10),
    legend.position = "right"  # you can also set "right", etc.
  )

# Apply the theme to all plots
plots_all <- list(
  s1_EC.S, s2_EC.S, s3_EC.S, s4_EC.S, s5_EC.S,
  s1_EC.p, s2_EC.p, s3_EC.p, s4_EC.p, s5_EC.p
)

plots_all <- lapply(plots_all, function(p) p + my_theme)

# Reassign them
s1_EC.S <- plots_all[[1]]
s2_EC.S <- plots_all[[2]]
s3_EC.S <- plots_all[[3]]
s4_EC.S <- plots_all[[4]]
s5_EC.S <- plots_all[[5]]
s1_EC.p <- plots_all[[6]]
s2_EC.p <- plots_all[[7]]
s3_EC.p <- plots_all[[8]]
s4_EC.p <- plots_all[[9]]
s5_EC.p <- plots_all[[10]]

# Top row: EC.S with scenario titles
top_row <- (
  (s1_EC.S + ggtitle(expression(atop("Scenario 1",paste(
    r[1], " = ", r[2], " = 1, ",
    s[1], " = ", s[2], " = 1"))))) |
    (s2_EC.S + ggtitle(expression(atop("Scenario 2",paste(
      r[1], " = ", r[2], " = 1, ",
      s[1], " = 1, ",
      s[2], " = 2"))))) |
    (s3_EC.S + ggtitle(expression(atop("Scenario 3",paste(
      r[1], " = ", r[2], " = 1, ",
      s[1], " = ", s[2], " = 2"))))) |
    (s4_EC.S + ggtitle(expression(atop("Scenario 4",paste(
      r[1], " = 1, ",
      r[2], " = 2, ",
      s[1], " = ", s[2], " = 1"))))) |
    (s5_EC.S + ggtitle(expression(atop("Scenario 5",paste(
      r[1], " = ", r[2], " = 2, ",
      s[1], " = ", s[2], " = 1")))))
)
# Bottom row: EC.P
bottom_row <- (
  s1_EC.p |
    s2_EC.p |
    s3_EC.p |
    s4_EC.p |
    s5_EC.p
)

# Combine rows and collect legends
final_plot <- (top_row / bottom_row) + 
  plot_layout(guides = "collect")
