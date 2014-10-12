plotArjas<-function(coxmod,data,time,event,censor="right"){

baseMod<-basehaz(coxmod,centered=TRUE)

survobj=Surv(data[,time],as.numeric(data[,event]),type=censor)
survmodObs=survfit(survobj~1, data=data)
lifetableObs=summary(survmodObs,time=baseMod$time)
cumObs<-cumsum(lifetableObs$n.event)
cumHaz=cumsum(baseMod$hazard)

#survmodExp=survfit(coxmod)
#lifetableExp=summary(survmodExp,time=baseMod$time)
#cumExp<-cumsum(lifetableExp$n.event)

maxLim=ceiling(cumObs[length(cumObs)]/10)*10
plot(stepfun(cumObs,c(0,cumHaz)),pty="sq",asp=1,cex=0.5,pch=16,xlim=c(0,maxLim),ylim=c(0,maxLim),
		main="Arjas plot",ylab="Expected number of events",xlab="Observed number of events")
abline(1,1)

}
