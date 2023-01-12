
#args=(commandArgs(TRUE))
#25 genetic variants are assumed to be candidate IVs and three scenarios
#are considered

# (1) Balanced pleiotropy - pleiotropic effects are equally likely to be positive
# as negative; (2) Directional pleiotropy - only positive pleiotropic effects are simulated; (3)
# Directional pleiotropy - pleiotropic effects are via a confounder. #InSIDE assumption is only valid in scenarios (1) and (2).


# The validity status of a genetic variant is determined by a random draw for each variant.
# where probability of being an invalid variant is varied in 3 scenarios with probabilities taken
# as 0.1, 0.2, and 0.3. Further, cases with 10,000 and 20,000 participants are considered.
# 10,000 simulated datasets are generated for each scenario in a two-sample setting with two
# values of the causal effect (β = 0, null causal effect; β = 0.1, positive causal effect)
#####functions####
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


#Data generation 

##change simulations setting as necessary
################one sample###################
###########################N=10000###################
####null causal effect B=0
#probability of being an invalid variant=0.1
# args=(commandArgs(TRUE))
iter=10000
IVW=matrix(data=NA, nrow=iter, ncol=2)
WM=matrix(data=NA, nrow=iter, ncol=2)
PWM=matrix(data=NA, nrow=iter, ncol=2)
EGGER=matrix(data=NA, nrow=iter, ncol=2)
set.seed(123)
for(k in 1:iter){

print(k)
  prob_IV=0.9
  B=0 # causal effect is NULL
  N=20000 # number of cases
  
  eps_u=rnorm(N,0,1)  # independent error terms for generating data
  eps_x=rnorm(N,0,1)
  eps_y=rnorm(N,0,1)
  
  genetic_variant=25
  #scenario1- balanced pleiotropy
  # phi=rep(0,genetic_variant) # INSIDE Assumption is satisfied
  # alpha=runif(genetic_variant,0,0.2) # balanced 
  
  phi=runif(genetic_variant,-0.2,0.2)
    #rep(0,genetic_variant) # INSIDE Assumption is satisfied
  alpha=runif(genetic_variant,0,0.2)
     # balanced 
  # G=data.frame() # generating Gij
  # for( j in 1:25){
  #   #assign(paste('G',j,sep = ''),sample(c(0,1,2),N,prob=c(0.49,0.42,0.09),replace = T))
  # G=rbind(G,sample(c(0,1,2),N,prob=c(0.49,0.42,0.09),replace = T))
  #   }
  #   
  G=matrix(c(rbinom(25*N,2,0.3)),nrow=25,ncol=20000)
  IV_indicator=rbinom(25,1,prob_IV) # deciding what are the valid instruments
  
  phi[IV_indicator==1]=0 #making phi and gamma of valid instruments 0.
  alpha[IV_indicator==1]=0
  #the genetic effects on gamma j are drawn from uniform distribution between 0.03 abd 0.1
  
  genetic_variant=25  # number of genetic variants
  gamma=runif(genetic_variant,0.03,0.1) # gamma_j
  
  effect1=unname(apply(t(G), 1, function(x) sum(x*phi)))
  effect2=unname(apply(t(G), 1, function(x) sum(x*gamma)))
  effect3=unname(apply(t(G), 1, function(x) sum(x*alpha)))
  
  
  
  confounder=effect1 +eps_u
  exposure=effect2 +confounder+eps_x
  outcome=effect3+B*exposure+confounder+eps_y
  
  data1=as.data.frame(t(G))
  colnames(data1)=c(paste('gene',1:25,sep = ''))
  data1$exposure=exposure
  data1$outcome=outcome
  
  library(broom)
  
  list_exposure <- lapply(colnames(data1)[1:25], function(i)
    tidy(lm(as.formula(paste("exposure ~", i)), data = data1)))
  summary_exposure=as.data.frame( t(sapply(list_exposure, "[",2, )))
  
  list_outcome <- lapply(colnames(data1)[1:25], function(i)
    tidy(lm(as.formula(paste("outcome ~", i)), data = data1)))
  summary_outcome=as.data.frame( t(sapply(list_outcome, "[",2, )))
  
  
  betaYG=unlist(summary_outcome$estimate)
  betaXG=unlist(summary_exposure$estimate)
  sebetaYG=unlist(summary_outcome$std.error)
  sebetaXG=unlist(summary_exposure$std.error)
  
  betaIV = betaYG/betaXG # ratio estimates
  
  
  # WM estimate
  penalty = pchisq(weights*(betaIV-betaIVW)^2, df=1, lower.tail=FALSE)
  pen.weights = weights*pmin(1, penalty*20) # penalized weights
  betaWM = weighted.median(betaIV, weights) # weighted median estimate
  sebetaWM = weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, weights)
  # standard error
  betaPWM = weighted.median(betaIV, pen.weights) # penalized weighted median estimate
  sebetaPWM = weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, pen.weights)
  # standard error
  
  #IVW estimate
  betaIVW = summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$coef[1,1]
  sebetaIVW = summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$coef[1,2]/
    min(summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$sigma, 1)
  
  
  
  betaYG=betaYG*sign(betaXG)
  betaXG=abs(betaXG)
  sebetaYG=unlist(summary_outcome$std.error)
  sebetaXG=unlist(summary_exposure$std.error)
  
  betaEGGER = summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$coef[2,1]
  sebetaEGGER = summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$coef[2,2]/
    min(summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$sigma, 1)
  
  
 IVW[k,]=c(betaIVW,sebetaIVW)
  WM[k,]=c(betaWM,sebetaWM)
  PWM[k,]=c(betaPWM,sebetaPWM)
  EGGER[k,] =c(betaEGGER,sebetaEGGER)
  
  # final_Data=rbind.data.frame(betaIVW,sebetaIVW,betaWM,sebetaWM,betaPWM,sebetaPWM,betaEGGER,sebetaEGGER)
  # row.names(final_Data)<-NULL
  # 
  # write.table(final_Data,file=paste0('out_',args), row.names=F, col.names=F)
 }

result=cbind(IVW,WM,PWM,EGGER)
colnames(result)=c("IVW","IVWse","WM","WMse","PWM","PWMse","EGGER","EGGERse")
# write.csv(IVW,"IVW_s2_1_1.csv")
# write.csv(WM,"WM_s2_1_1.csv")
# write.csv(PWM,"PWM_s2_1_1.csv")
# write.csv(EGGER,"EGGER_s2_1_1.csv")

write.csv(result,"result_s3_1_1_2n.csv")