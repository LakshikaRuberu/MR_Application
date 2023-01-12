library(readxl)
result <- read.csv("result_s3_2_3_2n.csv")
#colnames(result)=c("IVW","IVWse","WM","WMse","PWM","PWMse","EGGER","EGGERse")
#result=result[-(10001:10002),]

result$IVWlow=result$IVW-1.96*result$IVWse
result$IVWhigh=result$IVW+1.96*result$IVWse

result$WMlow=result$WM-1.96*result$WMse
result$WMhigh=result$WM+1.96*result$WMse

result$PWMlow=result$PWM-1.96*result$PWMse
result$PWMhigh=result$PWM+1.96*result$PWMse

result$EGGERlow=result$EGGER-1.96*result$EGGERse
result$EGGERhigh=result$EGGER+1.96*result$EGGERse

result$IndIVW=ifelse(result$IVWlow<=0 &result$IVWhigh>=0,0,1)
result$IndWM=ifelse(result$WMlow<=0 &result$WMhigh>=0,0,1)
result$IndPWM=ifelse(result$PWMlow<=0 &result$PWMhigh>=0,0,1)
result$IndEGGER=ifelse(result$EGGERlow<=0 &result$EGGERhigh>=0,0,1)

round(colMeans(result[,-c(1,10:17)]),3)
