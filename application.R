###Application

#reading the dataset
data=read.csv('dodata.csv',header = T)


#sample size for SNP exposure association
n1=188577
#sample size for SNP outcome association
n2=86995


#computing the SE of the association betas for the three exposures of interest and outcome
data$ldlse=abs(data$ldlcbeta)/qt(data$ldlcp2/2,188577-2,lower.tail = F)
data$hdlse=abs(data$hdlcbeta)/qt(data$hdlcp1/2,188577-2,lower.tail = F)
data$trise=abs(data$tgbeta)/qt(data$tgp1/2,188577-2,lower.tail = F)
data$cadse=abs(data$chdbeta)/qt(data$chdp1/2,86995-2,lower.tail = F)

summary(data$cadse)
data$cadse[is.na(data$cadse)==T]=50
data$trise[is.na(data$trise)==T]=50

#subset of SNPs associated with the exposures at a p value 1*10^-8
data_LDL1=subset(data,data$ldlcp2< 0.00000001)
data_HDL1=subset(data,data$hdlcp1< 0.00000001)
data_tri1=subset(data,data$tgp1< 0.00000001)


#subset of SNPs most associated with the exposure of interest
data_LDL1$minp=pmin(data_LDL1$ldlcp2,data_LDL1$hdlcp1,data_LDL1$tgp1)
data_LDL2=subset(data_LDL1,data_LDL1$ldlcp2== data_LDL1$minp)

data_HDL1$minp=pmin(data_HDL1$ldlcp2,data_HDL1$hdlcp1,data_HDL1$tgp1)
data_HDL2=subset(data_HDL1,data_HDL1$hdlcp1== data_HDL1$minp)

data_tri1$minp=pmin(data_tri1$ldlcp2,data_tri1$hdlcp1,data_tri1$tgp1)
data_tri2=subset(data_tri1,data_tri1$tgp1==data_tri1$minp)





#############################functions to perform MR methods###############

weighted.median <- function(betaIV.in, weights.in) {
  betaIV.order = betaIV.in[order(betaIV.in)]
  weights.order = weights.in[order(betaIV.in)]
  weights.sum = cumsum(weights.order)-0.5*weights.order
  weights.sum = weights.sum/sum(weights.order)
  below = max(which(weights.sum<0.5))
  weighted.est = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
    (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
  return(weighted.est) }

weighted.median.boot = function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in){
  med = NULL
  for(i in 1:1000){
    betaXG.boot = rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
    betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
    betaIV.boot = betaYG.boot/betaXG.boot
    med[i] = weighted.median(betaIV.boot, weights.in)
  }
  return(sd(med)) }


#plot(-data_LDL1$ldlcbeta,-data_LDL1$chdbeta)
#######################################################################
################HDL1################# #change the relevant dataset for other estimates.
betaYG=-data_HDL1$chdbeta
betaXG=-data_HDL1$hdlcbeta
sebetaYG=data_HDL1$cadse
sebetaXG=data_HDL1$hdlcp1
DF= length(betaYG)-1
betaIV = betaYG/betaXG # ratio estimates



weights = (sebetaYG/betaXG)^-2 # inverse-variance weights

#beta
betaIVW = summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$coef[1,1]
#SE
sebetaIVW = summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$coef[1,2]/
  min(summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$sigma, 1)
#p value
IVW_p = 2*(1-pt(abs(betaIVW/sebetaIVW),DF))


# WM methods estimate
penalty = pchisq(weights*(betaIV-betaIVW)^2, df=1, lower.tail=FALSE)
pen.weights = weights*pmin(1, penalty*20) # penalized weights
betaWM = weighted.median(betaIV, weights) # weighted median estimate
sebetaWM = weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, weights)
WM_p = 2*(1-pt(abs(betaWM/sebetaWM),DF))





betaPWM = weighted.median(betaIV, pen.weights) # penalized weighted median estimate
sebetaPWM = weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, pen.weights)
PWM_p = 2*(1-pt(abs(betaPWM/sebetaPWM),DF))







############################MR egger Regression####################
betaYG=data_HDL1$chdbeta*sign(data_HDL1$hdlcbeta) 
betaXG=abs(data_HDL1$hdlcbeta)
sebetaYG=data_HDL1$cadse
sebetaXG=data_HDL1$hdlcp1
DF= length(betaYG)-1
betaIV = betaYG/betaXG # ratio estimates







betaEGGER =  betaEGGER = summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$coef[2,1]
sebetaEGGER = summary(lm(betaYG~betaXG, weights=sebetaYG^-2))


betaEGGER
sebetaEGGER

