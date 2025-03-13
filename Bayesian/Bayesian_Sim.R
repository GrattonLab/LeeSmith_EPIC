#Loading necessary packages
install.packages('R2WinBUGS')
library('R2WinBUGS')
library(MASS)


#Set number of trials per subject per condition (TrialN)
#Set ratio of within subject variance to between subject variance (Var_Ratio)
#Set number of subjects (S)
TrialN<-50
Var_Ratio<-5
S<-100

#Simulating RT observation for conflict task (25 replications)
Sim_Data_List<-list()

STrialN<-S*TrialN 
WT<-(792*Var_Ratio)**.5
Sigma <- matrix(c(792,.075,.075,792),2,2)

w<-(mvrnorm(n = S, c(92.82, 92.82), Sigma, empirical = TRUE))

Sigma2 <- matrix(c(7000,.065,.065,7000),2,2)

m<-(mvrnorm(n = S, c(827.48,827.48), Sigma2, empirical = TRUE))


m1<-m[,1]
m2<-m[,2]

w1<-w[,1]
w2<-w[,2]

I1<-m1+(w1/2)
I2<-m2+(w2/2)

c1<-m1-(w1/2)
c2<-m2-(w2/2)

n1I<-matrix(NA,STrialN,1)
n1C<-matrix(NA,STrialN,1)
n2I<-matrix(NA,STrialN,1)
n2C<-matrix(NA,STrialN,1)

s1I<-rep(1:S,each=TrialN)
s1C<-rep(1:S,each=TrialN)
s2I<-rep(1:S,each=TrialN)
s2C<-rep(1:S,each=TrialN)

for (r in 1:25){ 
  for (i in 1:STrialN){
    mean1<-I1[s1I[i]]
    n1I[i,1]<-rnorm(1,mean1,WT)
  }
  
  for (i in 1:STrialN){
    mean2<-c1[s1C[i]]
    n1C[i,1]<-rnorm(1,mean2,WT)
  }
  
  for (i in 1:STrialN){
    mean3<-I2[s2I[i]]
    n2I[i,1]<-rnorm(1,mean3,WT)
  }
  
  
  for (i in 1:STrialN){
    mean4<-c2[s2C[i]]
    n2C[i,1]<-rnorm(1,mean4,WT)
  }
  
  n1I<-cbind.data.frame(n1I,s1I)
  n2I<-cbind.data.frame(n2I,s2I)
  n1C<-cbind.data.frame(n1C,s1C)
  n2C<-cbind.data.frame(n2C,s2C)
  
  CongCode<-rep(-.5,STrialN)
  IncCongCode<-rep(.5,STrialN)
  
  n1I$code<-IncCongCode
  n2I$code<-IncCongCode
  
  n1C$code<-CongCode
  n2C$code<-CongCode
  
  rt<-c(n1C$n1C,n1I$n1I,n2C$n2C,n2I$n2I)
  session<-rep(c(1,2),each=STrialN*2)
  sub<-c(n1C$s1C,n1I$s1I,n2C$s2C,n2I$s2I)
  cond<-c(n1C$code,n1I$code,n2C$code,n2I$code)
  Sim_Data_List[[r]]<-cbind.data.frame(rt,session,sub,cond)
}






#For each simulated dataset, use multilevel model to estimate:
#1)between subject variance
#2)within subject variance 
#Note that the code also calculates bias and precision for the variance estimates above
library(lme4)
MLM_sim_var<-matrix(NA,3,25)
MLM_sim_coef<-matrix(NA,2,25)
for (r in 1:25){
  AEL<-Sim_Data_List[[r]]
  AELs1<-AEL[AEL$session==1,]
  LG<-lmer(rt~cond+(1+cond|sub),data=AELs1,REML=FALSE)
  
  MLM_sim_var[,r] <- as.data.frame(VarCorr(LG))[c(1,2,4),4]
  MLM_sim_coef[,r]<-coef(summary(LG))[,1]
}

var_true<-c(7000,792,WT**2)
coef_true<-c(663,74)
MLM_varbias1<-matrix(NA,3,25)
MLM_varbias1[1,]<-MLM_sim_var[1,]-7000
MLM_varbias1[2,]<-MLM_sim_var[2,]-792
MLM_varbias1[3,]<-MLM_sim_var[3,]-WT**2

MLM_coefbias1<-matrix(NA,2,25)
MLM_coefbias1[1,]<-MLM_sim_coef[1,]-827.48
MLM_coefbias1[2,]<-MLM_sim_coef[2,]-92.82

MLM_varbias2<-apply(MLM_varbias1,1,mean)
MLM_coefbias2<-apply(MLM_coefbias1,1,mean)

MLM_varMAD1<-abs(MLM_varbias1)
MLM_varMAD2<-apply(MLM_varMAD1,1,mean)

MLM_coefMAD1<-abs(MLM_coefbias1)
MLM_coefMAD2<-apply(MLM_coefMAD1,1,mean)

MLM_bias_all<-c(MLM_coefbias2,MLM_varbias2)
MLM_MAD_all<-c(MLM_coefMAD2,MLM_varMAD2)

write.csv(MLM_bias_all,'MLM_bias_all.csv')
write.csv(MLM_MAD_all,'MLM_MAD_all.csv')

#Obtaining Bayesian estimates of subject-level congruency effects
#Note that for the Bayesian estimates, the within subject variance and between subject variance from the MLM are used as priors
Bayes_Est<-matrix(NA,100,25)


for (r in 1:25){
  z<-Sim_Data_List[[r]]
  y<-z[z$session==1,]
  mu_b0=MLM_sim_coef[1,r]
  mu_b1=MLM_sim_coef[2,r]
  prec_b0<-1/MLM_sim_var[1,r]
  prec_b1=1/MLM_sim_var[2,r]
  prec_error=1/MLM_sim_var[3,r]
  All_List<-list(y=as.matrix(y),N=100*TrialN*2,P=100,mu_b0=mu_b0,mu_b1=mu_b1,prec_b0=prec_b0,prec_b1=prec_b1,prec_error=prec_error)
  
  BayesResult <- bugs(data=All_List, parameters.to.save=c('b1','var_b1'), model.file='Cong.txt',
                    n.chains=3, n.iter=20000,n.burnin=10000,inits=NULL,
                    bugs.directory="C:/Users/Skip/Desktop/WB/winbugs14_full_patched/WinBUGS14")$summary
  Bayes_Est[,r]<-BayesResult[1:100,1]
}





#Calculating bias and precision values for Bayesian estimates of subject-level congruency effects:

wTrue<-c(w1)
bias1T<-Bayes_Est-wTrue

bias2T<-apply(bias1T,1,sum)
bias3T<-bias2T/25



mad1<-abs(bias1T)
mad2<-apply(mad1,1,mean)


write.csv(bias3T,'BAYES_bias.csv')
write.csv(mad2,'BAYES_mad.csv')


#Calculating subject-level congruency effects with a two-step frequentist approach:
Trad_Est<-matrix(NA,100,25)
for (r in 1:25){
  Everything<-(Sim_Data_List[[1]])
  
  
  
  single_reg_coefs1<-matrix(NA,100,2)
  single_reg_residvars1<-matrix(NA,100,1)
  for (p in 1:100){
    trim<-Everything[Everything$sub==p&Everything$session==1,]
    trim_reg<-lm(trim$rt~trim$cond)
    single_reg_coefs1[p,]<-trim_reg$coefficients
    single_reg_residvars1[p,1]<-var(trim_reg$residuals)
  }

  Trad_Est[,r]<-single_reg_coefs1[,2]
}

#Calculating bias and precision for the two-step frequentist estimates of congruency effect:
bias1T<-Trad_Est-w1

bias2T<-apply(bias1T,1,sum)
bias3T<-bias2T/25

mad1<-abs(bias1T)
mad2<-apply(mad1,1,mean)

var_cong<-apply(Trad_Est,2,var)
bias1T_var<-var_cong-5.3

bias2T_var<-sum(bias1T_var)
bias3T_var<-bias2T_var/25

mad_var<-mean(abs(bias1T_var))

bias_all<-c(bias3T,bias3T_var)
mad_all<-c(mad2,mad_var)

write.csv(bias_all,'Trad_bias.csv')
write.csv(mad_all,'Trad_mad.csv')






