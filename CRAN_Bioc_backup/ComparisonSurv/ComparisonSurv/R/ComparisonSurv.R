
##################################################################
######      Descriptive statistic for time-to-event data    ######
##################################################################
Descrip.method<-function(time,status,group,tau="observed",alpha=0.05){
  index0<-which(group==0)
  index1<-which(group==1)
  group0.samplesize<-length(index0)
  group1.samplesize<-length(index1)
  group0.event.size<-sum(status[index0])
  group1.event.size<-sum(status[index1])
  group0.censore.rate<-1-sum(status[index0])/group0.samplesize
  group1.censore.rate<-1-sum(status[index1])/group1.samplesize
  group0.max.event.time<-max(time[index0][which(status[index0]==1)])
  group1.max.event.time<-max(time[index1][which(status[index1]==1)])
  group0.max.time<-max(time[index0])
  group1.max.time<-max(time[index1])

  kmfit<-survfit(Surv(time,status==1)~group)
  s<-summary(kmfit,alpha=alpha)
  su<-s$table
  group0.mean.time<-su[1,5]
  group0.mean.se<-su[1,6]
  group0.mean.lower95<-group0.mean.time-qnorm(1-alpha/2)*sqrt(su[1,6])
  group0.mean.upper95<-group0.mean.time+qnorm(1-alpha/2)*sqrt(su[1,6])
  group1.mean.time<-su[2,5]
  group1.mean.se<-su[2,6]
  group1.mean.lower95<-group1.mean.time-qnorm(1-alpha/2)*sqrt(su[2,6])
  group1.mean.upper95<-group1.mean.time+qnorm(1-alpha/2)*sqrt(su[2,6])

  sq<-quantile(kmfit)
  group0.median.time<-su[1,7]
  group0.median.lower95<-su[1,8]
  group0.median.upper95<-su[1,9]
  group1.median.time<-su[2,7]
  group1.median.lower95<-su[2,8]
  group1.median.upper95<-su[2,9]

  group0.time.25<-sq$quantile[,"25"][1]
  group0.time.25.lower95<- sq$lower[,"25"][1]
  group0.time.25.upper95<- sq$upper[,"25"][1]
  group1.time.25<-sq$quantile[,"25"][2]
  group1.time.25.lower95<- sq$lower[,"25"][2]
  group1.time.25.upper95<- sq$upper[,"25"][2]

  group0.time.75<-sq$quantile[,"75"][1]
  group0.time.75.lower95<- sq$lower[,"75"][1]
  group0.time.75.upper95<- sq$upper[,"75"][1]
  group1.time.75<-sq$quantile[,"75"][2]
  group1.time.75.lower95<- sq$lower[,"75"][2]
  group1.time.75.upper95<- sq$upper[,"75"][2]

  strata<-s$strata

  if(tau=="observed"){
    tau=min(group0.max.time,group1.max.time)
  }

  if(tau=="event"){
    tau=min(group0.max.event.time,group1.max.event.time)
  }

  tau1<-tau
  if(tau1>max(group0.max.time,group1.max.time)) stop("tau is larger than the maximum of observed time in both of the two groups")

  rm<-rmst2(time,status,arm=group,tau1,alpha=alpha)
  #  plot(rm)
  group0.RMST<-rm$RMST.arm0$rmst[[1]]
  group0.RMST.se<-rm$RMST.arm0$rmst[[2]]
  group0.RMST.lower95<-rm$RMST.arm0$rmst[[3]]
  group0.RMST.upper95<-rm$RMST.arm0$rmst[[4]]
  group1.RMST<-rm$RMST.arm1$rmst[[1]]
  group1.RMST.se<-rm$RMST.arm1$rmst[[2]]
  group1.RMST.lower95<-rm$RMST.arm1$rmst[[3]]
  group1.RMST.upper95<-rm$RMST.arm1$rmst[[4]]

  Y<-list()
  Y$result.summary=data.frame(sample.size=c(group0.samplesize,group1.samplesize),
                              event.num=c(group0.event.size,group1.event.size),
                              censoring.rate=c(group0.censore.rate,group1.censore.rate),
                              max.observed.time=c(group0.max.time,group1.max.time),
                              max.event.time=c(group0.max.event.time,group1.max.event.time))
  rownames(Y$result.summary)<-c("group=0","group=1")

  Y$result.mean=data.frame(   mean.time=c(group0.mean.time,group1.mean.time),
                              se=c(group0.mean.se,group1.mean.se),
                              mean.lower95=c(group0.mean.lower95,group1.mean.lower95),
                              mean.upper95=c(group0.mean.upper95,group1.mean.upper95))
  rownames(Y$result.mean)<-c("group=0","group=1")
  colnames(Y$result.mean)<-c("Est.","se",paste("lower .", round((1 -  alpha) * 100, digits = 0), sep = ""),
                             paste("upper .",round((1 - alpha) * 100, digits = 0), sep = ""))

  Y$result.quantile=data.frame(time.25=c(group0.time.25,group1.time.25),
                               time.25.lower95=c(group0.time.25.lower95,group1.time.25.lower95),
                               time.25.upper95=c(group0.time.25.upper95,group1.time.25.upper95),
                               median.time=c(group0.median.time,group1.median.time),
                               median.lower95=c(group0.median.lower95,group1.median.lower95),
                               median.upper95=c(group0.median.upper95,group1.median.upper95),
                               time.75.time=c(group0.time.75,group1.time.75),
                               time.75.lower95=c(group0.time.75.lower95,group1.time.75.lower95),
                               time.75.upper95=c(group0.time.75.upper95,group1.time.75.upper95))
  rownames(Y$result.quantile)<-c("group=0","group=1")
  colnames(Y$result.quantile)<-c("Est.25",paste("lower .", round((1 -  alpha) * 100, digits = 0), sep = ""),
                                 paste("upper .",round((1 - alpha) * 100, digits = 0), sep = ""),
                                 "Est.50",paste("lower .", round((1 -  alpha) * 100, digits = 0), sep = ""),
                                 paste("upper .",round((1 - alpha) * 100, digits = 0), sep = ""),
                                 "Est.75",paste("lower .", round((1 -  alpha) * 100, digits = 0), sep = ""),
                                 paste("upper .",round((1 - alpha) * 100, digits = 0), sep = ""))

  Y$tau<-tau1
  Y$result.RMST=data.frame(   RMST=c(group0.RMST,group1.RMST),
                              RMST.se=c(group0.RMST.se,group1.RMST.se),
                              RMST.lower95=c(group0.RMST.lower95,group1.RMST.lower95),
                              RMST.upper95=c(group0.RMST.upper95,group1.RMST.upper95))
  rownames(Y$result.RMST)<-c("group=0","group=1")
  colnames(Y$result.RMST)<-c("Est.","se",paste("lower .", round((1 -  alpha) * 100, digits = 0), sep = ""),
                             paste("upper .",round((1 - alpha) * 100, digits = 0), sep = ""))

  Y
}

Descriptive.stat<-function(time,status,group,tau="observed",alpha=0.05){

  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))!=2) stop("There must be two groups")
  if (all(group%in%c(0,1))==FALSE)  stop("Group must be 0 or 1")

  Y<-Descrip.method(time,status,group,tau=tau,alpha=alpha)
  tau<-Y$tau
  print(Y)
}

######################################################################
######                 Test overall homogeneity                 ######
######################################################################
Bre_TW<-function(time,status,group,weight=1){

  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))!=2)     stop("There must be two groups")
  if (all(group%in%c(1,2))==FALSE)  stop("Group must be 1 or 2")
  if((weight!=1)&(weight!=0.5))     stop("Weight is an invalid value")

  fit<-survfit(Surv(time,status)~group)
  t<-fit$time[1:fit$strata[1]]
  t2<-fit$time[(fit$strata[1]+1):(fit$strata[1]+fit$strata[2])]
  n1<-fit$n.risk[1:fit$strata[1]]
  n2<-fit$n.risk[(fit$strata[1]+1):(fit$strata[1]+fit$strata[2])]
  d1<-fit$n.event[1:fit$strata[1]]
  d2<-fit$n.event[(fit$strata[1]+1):(fit$strata[1]+fit$strata[2])]

  df1<-data.frame(t,n1,d1)
  df2<-data.frame(t2,n2,d2)
  df<-merge(df1,df2,all=TRUE,by.x="t",by.y="t2")

  df$d1[is.na(df$d1)]<-0
  df$d2[is.na(df$d2)]<-0

  new<-df[which((df$d1!=0)|(df$d2!=0)),]
  tn1<-new$t[is.na(new$n1)]
  tn2<-new$t[is.na(new$n2)]

  num1<-0
  num2<-0
  for(i in 1:length(tn1))
  {num1[i]<-sum((time[group==1])>=tn1[i])}
  for(j in 1:length(tn2))
  {num2[j]<-sum((time[group==2])>=tn2[j])}
  new[which(is.na(new$n1)),2]<-num1
  new[which(is.na(new$n2)),4]<-num2

  risk<-new$n1+new$n2
  event<-new$d1+new$d2
  e1<-new$n1*event/risk
  e2<-new$n2*event/risk
  diff<-new$d1-e1
  V<-(new$n1/risk)*(1-new$n1/risk)*((risk-event)/(risk-1))*event
  all<-data.frame(t=new$t,n1=new$n1,d1=new$d1,e1,n2=new$n2,d2=new$d2,e2,risk,event,diff,V)

  V[is.nan(V)]<-0
  if(weight==1)
  {
    Stat1<-(sum(risk*diff))/sqrt(sum((risk^2)*V))
    pvalue1<-2*(1-pnorm(abs(Stat1)))
    result1<-data.frame(Method="Gehan-Wilcoxon",Weight="Yi",Statistic=Stat1,pvalue=pvalue1)
    return(result1)
  }
  if(weight==0.5)
  {
    Stat2<-(sum(sqrt(risk)*diff))/sqrt(sum(risk*V))
    pvalue2<-2*(1-pnorm(abs(Stat2)))
    result2<-data.frame(Method="Tarone-Ware",Weight="sqrt(Yi)",Statistic=Stat2,pvalue=pvalue2)
    return(result2)
  }
}

WKM<-function(time,status,group){
  dataset<-data.frame(time,status,group+1)
  fitall<-survfit(Surv(time,status)~1)
  df1<-survfit(Surv(time,status)~1,data=dataset[dataset$group==1,])
  df2<-survfit(Surv(time,status)~1,data=dataset[dataset$group==2,])
  n1<-df1$n
  n2<-df2$n
  cf1<-survfit(Surv(time,1-status)~1,data=dataset[dataset$group==1,])
  cf2<-survfit(Surv(time,1-status)~1,data=dataset[dataset$group==2,])
  findt<-function(t){
    if(t>max(df1$time))
    {s1t<-summary(df1,times=max(df1$time))$surv
    c1t<-summary(cf1,times=max(cf1$time))$surv
    s2t<-summary(df2,times=t)$surv
    c2t<-summary(cf2,times=t)$surv
    }
    else if(t>max(df2$time))
    {s1t<-summary(df1,times=t)$surv
    c1t<-summary(cf1,times=t)$surv
    s2t<-summary(df2,times=max(df2$time))$surv
    c2t<-summary(cf2,times=max(cf2$time))$surv
    }
    else
    {s1t<-summary(df1,times=t)$surv
    s2t<-summary(df2,times=t)$surv
    c1t<-summary(cf1,times=t)$surv
    c2t<-summary(cf2,times=t)$surv
    }
    data.frame(t,s=summary(fitall,times=t)$surv,s1t,s2t,c1t,c2t)
  }
  dat<-numeric(0)
  for(i in fitall$time)
  {dat<-rbind(dat,findt(i))}
  dat<-subset(dat,((dat$c1t>0)&(dat$c2t>0)&(dat$s1t>0)&(dat$s2t>0)))
  s_=c(1,dat$s[-length(dat$s)])
  c1t_=c(1,dat$c1t[-length(dat$c1t)])
  c2t_=c(1,dat$c2t[-length(dat$c2t)])
  dt=c(diff(dat$t,1),0) #dt=ti-t(i-1)
  ds=dat$s-s_
  w=((n1+n2)*c1t_*c2t_)/((n1*c1t_)+(n2*c2t_))
  dat<-cbind(dat,c1t_,c2t_,dt,w,s_,ds)
  delmax<-dat[1:(length(dat$t)-1),]
  WKM<-with(delmax,sqrt(n1*n2/(n1+n2))*sum((w*(s2t-s1t)*dt)))
  mult<-delmax$w*delmax$s*delmax$dt
  inter<-(cumsum(mult[length(mult):1]))[length(mult):1]
  var22=-sum(with(delmax,(inter^2)*ds/(s*s_*w)))
  Statistic<-WKM/sqrt(var22)
  Pvalue<-2*(1-pnorm(abs(Statistic)))
  result<-data.frame(Method="Weighted Kaplan-Meier",Statistic,Pvalue)
}

LWmethod<-function(time,status,group){
  group<-group+1
  fit<-survfit(Surv(time,status)~group)
  t<-fit$time[1:fit$strata[1]]
  t2<-fit$time[(fit$strata[1]+1):(fit$strata[1]+fit$strata[2])]
  n1<-fit$n.risk[1:fit$strata[1]]
  n2<-fit$n.risk[(fit$strata[1]+1):(fit$strata[1]+fit$strata[2])]
  d1<-fit$n.event[1:fit$strata[1]]
  d2<-fit$n.event[(fit$strata[1]+1):(fit$strata[1]+fit$strata[2])]
  df1<-data.frame(t,n1,d1)
  df2<-data.frame(t2,n2,d2)
  df<-merge(df1,df2,all=T,by.x="t",by.y="t2")
  df$d1[is.na(df$d1)]<-0
  df$d2[is.na(df$d2)]<-0
  new<-df[which((df$d1!=0)|(df$d2!=0)),]
  tn1<-new$t[is.na(new$n1)]
  tn2<-new$t[is.na(new$n2)]
  num1<-0
  num2<-0
  for(i in 1:length(tn1))
  {num1[i]<-sum((time[group==1])>=tn1[i])}
  for(j in 1:length(tn2))
  {num2[j]<-sum((time[group==2])>=tn2[j])}
  new[which(is.na(new$n1)),2]<-num1
  new[which(is.na(new$n2)),4]<-num2
  risk<-new$n1+new$n2
  event<-new$d1+new$d2
  e1<-new$n1*event/risk
  e2<-new$n2*event/risk
  all<-data.frame(t=new$t,n1=new$n1,d1=new$d1,e1,n2=new$n2,d2=new$d2,e2,risk,event)
  delta<-sum((all$d1-e1)^2)
  temp1<-(all$n1*all$n2*event*(risk-event))/((risk^2)*(risk-1))
  temp1[risk<=1]<-0
  Edelta<-sum(temp1)
  V1<-all$n1*all$n2*event*(risk-event)/((risk^2)*(risk-1))
  V1[risk<=1]<-0
  E2<-V1+e1^2
  fac3<-event*(event-1)*(event-2)*all$n1*(all$n1-1)*(all$n1-2)/(risk*(risk-1)*(risk-2)*factorial(3))
  fac3[risk<=2]<-0
  E3<-3*E2-2*e1+factorial(3)*fac3
  fac4<-event*(event-1)*(event-2)*(event-3)*all$n1*(all$n1-1)*(all$n1-2)*(all$n1-3)/(risk*(risk-1)*(risk-2)*(risk-3)*factorial(4))
  fac4[risk<=3]<-0
  E4<-6*E3-11*E2+6*e1+factorial(4)*fac4
  Vdelta<-sum(E4-4*E3*e1+6*E2*e1^2-3*e1^4-V1^2)
  Statistic<-(delta-Edelta)/sqrt(Vdelta)
  Pvalue<-2*(1-pnorm(abs(Statistic)))
  result<-data.frame(Method="Linwang",Edelta=Edelta,Vdelta=Vdelta,Statistic,Pvalue)
  result
}

Linxumethod<-function(time,status,group,side=c("one.sided","two.sided"),rou=0.5){
  group<-group+1
  if ((rou<0)|(rou>=1)) stop("rou must between 0 and 1")
  dataset<-data.frame(time,status,group)
  fitall<-survfit(Surv(time,status)~1)
  df1<-survfit(Surv(time,status)~1,data=dataset[dataset$group==1,])
  df2<-survfit(Surv(time,status)~1,data=dataset[dataset$group==2,])
  n1<-df1$n
  n2<-df2$n
  tds1<-data.frame(t=df1$time,d1=df1$n.event,s1t=df1$surv)
  tds2<-data.frame(t2=df2$time,d2=df2$n.event,s2t=df2$surv)
  tds<-merge(tds1,tds2,all=T,by.x="t",by.y="t2")
  tds$d1[is.na(tds$d1)]<-0
  tds$d2[is.na(tds$d2)]<-0
  for(i in 2:nrow(tds)){
    if(is.na(tds$s1t[i])) tds$s1t[i]<-tds$s1t[i-1]
    if(is.na(tds$s2t[i])) tds$s2t[i]<-tds$s2t[i-1]
  }
  tds$s1t[is.na(tds$s1t)]<-1
  tds$s2t[is.na(tds$s2t)]<-1
  lastevent<-c((tds1$d1[nrow(tds1)]),(tds2$d2[nrow(tds2)]))
  lasttime<- c((tds1$t[nrow(tds1)]),(tds2$t2[nrow(tds2)]))
  if(all(lastevent==0)) tau=min(lasttime)
  if(any(lastevent==0)&any(lastevent!=0))
    tau=max((lasttime[1]*(1-lastevent[1])),(lasttime[2]*(1-lastevent[2])))
  if(all(lastevent!=0)) tau=max(lasttime)
  table1<-subset(tds,(tds$d1!=0)|(tds$d2!=0))
  table1<-rbind(table1,tds[tds$t==tau,])
  dt<-c(diff(table1$t,1),0) #calculate the t(i+1)-ti
  delta<-with(table1,sum(abs(s1t-s2t)*dt)) #delta=5.1201
  vs1t<-with(df1,(surv^2)*cumsum(n.event/(n.risk*(n.risk-n.event))))
  vs2t<-with(df2,(surv^2)*cumsum(n.event/(n.risk*(n.risk-n.event))))
  vs1t[is.nan(vs1t)]<-0
  vs2t[is.nan(vs2t)]<-0
  vs1table<-data.frame(time1=df1$time,vs1t)
  vs2table<-data.frame(time2=df2$time,vs2t)
  vstable<-merge(vs1table,vs2table,all=TRUE,by.x="time1",by.y="time2")
  for(k in 2:nrow(vstable))
  {
    if(is.na(vstable$vs1t[k])) vstable$vs1t[k]<-vstable$vs1t[k-1]
    if(is.na(vstable$vs2t[k])) vstable$vs2t[k]<-vstable$vs2t[k-1]
  }
  vstable$vs1t[is.na(vstable$vs1t)]<-0
  vstable$vs2t[is.na(vstable$vs2t)]<-0
  table2<-merge(table1,vstable,by.x="t",by.y="time1")
  Edelta<-sum(sqrt((2/pi)*(table2$vs1t+table2$vs2t))*dt)
  Vardelta1<-sum((1-2/pi)*(table2$vs1t+table2$vs2t)*(dt^2))
  Vardelta2<-0
  for(m in 1:(nrow(table2)-2)){
    for(n in (m+1):(nrow(table2)-1)){
      Vardelta2=Vardelta2+with(table2,2*rou*(t[m+1]-t[m])*(t[n+1]-t[n])*(1-2/pi)*sqrt((vs1t[m]+vs2t[m])*(vs1t[n]+vs2t[n])))
    }
  }
  deltastar<-(delta-Edelta)/sqrt(Vardelta1+Vardelta2)
  if(side=="one.sided") Pvalue<-1-pnorm(deltastar)
  if(side=="two.sided") Pvalue<-2*(1-pnorm(abs(deltastar)))
  result<-data.frame(Method="LinXu method",delta,Edelta,Vardelta=Vardelta1+Vardelta2,
                     side,rou,Statistic=deltastar,Pvalue)
  result
}

Overall.test<-function(time,status,group,tau=NULL,nperm=500,seed=12345){
  set.seed(seed)

  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))!=2) stop("There must be two groups")
  if (all(group%in%c(0,1))==FALSE)  stop("Group must be 0 or 1")

  #PH
  fit <- coxph(Surv(time,status==1)~group)
  temp<- cox.zph(fit)
  ph_s<-temp$table["GLOBAL","chisq"]
  ph_p<-temp$table["GLOBAL","p"]

  #twostage
  twostage<-twostage(time,status,group,nperm) #twostage sample size 1000
  twostage_p<-twostage[[3]]
  twostage_s<-"NA"

  #logrank
  logrank<-survdiff(Surv(time,status)~group,rho=0)
  logrank_s<-logrank$chisq
  logrank_p<-1-pchisq(logrank_s,1)

  #Breslow-TW test
  GW<-Bre_TW(time,status,group+1,weight=1)
  GW_s<-GW$Statistic
  GW_p<-GW$pvalue

  TW<-Bre_TW(time,status,group+1,weight=0.5)
  TW_s<-TW$Statistic
  TW_p<-TW$pvalue

  #WKM
  a<-WKM(time,status,group)
  WKM_s<-a$Statistic
  WKM_p<-a$Pvalue

  #Linwang Method
  b<-LWmethod(time,status,group)
  LWmethod_s<-b$Statistic
  LWmethod_p<-b$Pvalue

  #lin-xu Method
  #rou: the correlation coefficient
  ccc<-Linxumethod(time,status,group,side="two.sided",rou=0.5)
  Linxu_s<-ccc$Statistic
  Linxu_p<-ccc$Pvalue

  #lin-xu permutation
  size1<-sum(group==0)
  size2<-sum(group==1)
  sta<-c()
  st2<-c()

  for(l in 1:nperm){
    group1<-sample(c(rep(0,size1),rep(1,size2)))
    sta[l]<-Linxumethod(time,status,group1,side="two.sided",rou=0.5)$Statistic
  }
  dd<-Linxumethod(time,status,group,side="two.sided",rou=0.5)$Statistic
  st2.PV<-sum(as.numeric(abs(sta)>abs(dd)))/nperm
  linxu_p2_s<-Linxumethod(time,status,group,side="two.sided",rou=0.5)$Statistic
  linxu_p2_p<-st2.PV

  #RMST
  rm<-rmst2(time,status,arm=group,tau)
  rmst_s_tau_15<-rm$unadjusted.result[1,1]
  rmst_p<-rm$unadjusted.result[1,4]
  tau<-rm$tau

  message("\n PH assumption yield a result of: Pvalue =", round(ph_p,5), ", chisq =", round(ph_s,5),".\n")

  result<-data.frame(method=c("Log-rank","Gehan-Wilcoxon","Tarone-Ware","Weighted KM"
                              ,"ABS","ABS permutation","Two-stage","Squared differences",paste("RMST (tau=",tau,")")),
                     statistic=c(round(logrank_s,5),round(GW_s,5),round(TW_s,5),round(WKM_s,5)
                                 ,round(Linxu_s,5),round(linxu_p2_s,5),twostage_s,round(LWmethod_s,5),round(rmst_s_tau_15,5)),
                     pvalue=c(round(logrank_p,5),round(GW_p,5),round(TW_p,5),round(WKM_p,5)
                              ,round(Linxu_p,5),round(linxu_p2_p,5),round(twostage_p,5),round(LWmethod_p,5),round(rmst_p,5)))
  print(result)
}



######################################################################
######           Test homogeneity at a fixed time point         ######
######################################################################

Fixpoint<-function(time,status,group,t0){
  kmfit1<-survfit(Surv(time,status==1)~1)
  summ<-summary(kmfit1,time=c(t0))
  surv<-summ$surv
  Vs<-summ$std.err
  std<-Vs^2/surv^2
  kmfit<-survfit(Surv(time,status==1)~group)
  summ1<-summary(kmfit,time=c(t0))
  surv1<-summ1$surv[1]
  Vs1<-summ1$std.err[1]
  std1<-(Vs1)^2/surv1^2
  surv2<-summ1$surv[2]
  Vs2<-summ1$std.err[2]
  std2<-(Vs2)^2/surv2^2
  n1<-summ1$n[1]
  n2<-summ1$n[2]
  n<-summ$n

  X1<-((surv1-surv2)^2)/(surv1^2*std1+surv2^2*std2)
  pX1<-1-pchisq(X1,df=1)
  ci11<-surv1-1.96*sqrt(surv1^2*std1)
  ci12<-surv1+1.96*sqrt(surv1^2*std1)
  ci13<-surv2-1.96*sqrt(surv2^2*std2)
  ci14<-surv2+1.96*sqrt(surv2^2*std2)

  X2<-((log(surv1)-log(surv2))^2)/(std1+std2)
  pX2<-1-pchisq(X2,df=1)
  ci21<-surv1^(1/(exp(1.96*std1/log(surv1))))
  ci22<-surv1^(exp(1.96*std1/log(surv1)))
  ci23<-surv2^(1/(exp(1.96*std2/log(surv2))))
  ci24<-surv2^(exp(1.96*std2/log(surv2)))

  X3<-((log(-log(surv1))-log(-log(surv2)))^2)/(std1/(log(surv1)^2)+std2/(log(surv2)^2))
  pX3<-1-pchisq(X3,df=1)
  ci31<-exp(-exp((log(-log(surv1)))-1.96*std1*log(surv1)))
  ci32<-exp(-exp((log(-log(surv1)))+1.96*std1*log(surv1)))
  ci33<-exp(-exp((log(-log(surv2)))-1.96*std2*log(surv2)))
  ci34<-exp(-exp((log(-log(surv2)))+1.96*std2*log(surv2)))

  v1<-surv1*std1/(4*(1-surv1))
  v2<-surv2*std2/(4*(1-surv2))
  X4<-((asin(sqrt(surv1))-asin(sqrt(surv2)))^2)/(v1+v2)
  pX4<-1-pchisq(X4,df=1)
  ci41<-sin(max(0,asin(surv1^(1/2))-0.5*1.96*std1*((surv1/(1-surv1))^(1/2))))^2
  ci42<-sin(min(pi/2,asin(surv1^(1/2))+0.5*1.96*std1*((surv1/(1-surv1))^(1/2))))^2
  ci43<-sin(max(0,asin(surv2^(1/2))-0.5*1.96*std2*((surv2/(1-surv2))^(1/2))))^2
  ci44<-sin(min(pi/2,asin(surv2^(1/2))+0.5*1.96*std2*((surv2/(1-surv2))^(1/2))))^2

  X5<-((log(surv1/(1-surv1))-log(surv2/(1-surv2)))^2)/(std1/((1-surv1)^2)+std2/((1-surv2)^2))
  pX5<-1-pchisq(X5,df=1)
  L1<-exp(1.96*std1/(1-surv1))
  ci51<-surv1/((L1-1)*surv1+L1)
  ci52<-L1*surv1/((L1-1)*surv1+1)
  L2<-exp(1.96*std2/(1-surv2))
  ci53<-surv2/((L2-1)*surv2+L2)
  ci54<-L2*surv2/((L2-1)*surv2+1)
  result<-list()
  result$est.g0<-data.frame(method=c("Naive","Log","Cloglog","Arcsin-square","Logit"),
                            t0=rep(t0,5),
                            est=rep(surv1,5),
                            "lower 95 CI"=c(ci11,ci21,ci31,ci41,ci51),
                            "upper 95 CI"=c(ci12,ci22,ci32,ci42,ci52))
  result$est.g1<-data.frame(method=c("Naive","Log","Cloglog","Arcsin-square","Logit"),
                            t0=rep(t0,5),
                            est=rep(surv2,5),
                            "lower 95 CI"=c(ci13,ci23,ci33,ci43,ci53),
                            "upper 95 CI"=c(ci14,ci24,ci34,ci44,ci54))
  result$test<-data.frame(method=c("Naive","Log","Cloglog","Arcsin-square","Logit"),
                          statistic=c(round(X1,5),round(X2,5),round(X3,5),round(X4,5),round(X5,5)),
                          pvalue=c(round(pX1,5),round(pX2,5),round(pX3,5),round(pX4,5),round(pX5,5)))
  print(result)
}

Fixpoint.test<-function(time,status,group,t0){

  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))<2) stop("There must be two or more groups")
  if(is.na(t0)==TRUE) stop("t0 is null,please select a time point")
  if(t0<=min(time[which(status==1)])) stop("t0 is too small,please select a time point more than the minimum non-censored timepoint")
  if(t0>=max(time[which(status==1)])) stop("t0 is too large,please select a time point less than the maximum non-censored timepoint")
  if (all(group%in%c(0,1))==FALSE)  stop("Group must be 0 or 1")

  Fixpoint(time,status,group,t0)
}



######################################################################
######           Test homogeneity before a time point           ######
######################################################################

Short<-function(time,status,group,t0){
  index<-which(time>t0)
  group<-group+1
  time[index]<-t0
  status[index]<-0
  dataset<-data.frame(time,status,group)
  t1<-dataset$time[which(dataset$status==1 & dataset$group==1)]
  t2<-dataset$time[which(dataset$status==1 & dataset$group==2)]
  if(t0<max(min(t1),min(t2)))
  {print("t0 is too small,please select another time greater than both minimum time")}

  logrank<-survdiff(Surv(time,status)~group)
  logrank.stat<-logrank$chisq
  logrank.P<-1-pchisq(logrank.stat,1)

  fitall<-survfit(Surv(time,status)~1)
  df1<-survfit(Surv(time,status)~1,data=dataset[dataset$group==1,])
  df2<-survfit(Surv(time,status)~1,data=dataset[dataset$group==2,])
  n1<-df1$n
  n2<-df2$n
  cf1<-survfit(Surv(time,1-status)~1,data=dataset[dataset$group==1,])
  cf2<-survfit(Surv(time,1-status)~1,data=dataset[dataset$group==2,])
  findt<-function(t){
    if(t>max(df1$time)){
      s1t<-summary(df1,times=max(df1$time))$surv
      c1t<-summary(cf1,times=max(cf1$time))$surv
      s2t<-summary(df2,times=t)$surv
      c2t<-summary(cf2,times=t)$surv}
    else if(t>max(df2$time)){
      s1t<-summary(df1,times=t)$surv
      c1t<-summary(cf1,times=t)$surv
      s2t<-summary(df2,times=max(df2$time))$surv
      c2t<-summary(cf2,times=max(cf2$time))$surv}
    else{
      s1t<-summary(df1,times=t)$surv
      s2t<-summary(df2,times=t)$surv
      c1t<-summary(cf1,times=t)$surv
      c2t<-summary(cf2,times=t)$surv}
    data.frame(t,s=summary(fitall,times=t)$surv,s1t,s2t,c1t,c2t)
  }
  dat<-numeric(0)
  for(i in fitall$time){
    dat<-rbind(dat,findt(i))}
  dat<-subset(dat,((dat$c1t>0)&(dat$c2t>0)&(dat$s1t>0)&(dat$s2t>0)))
  s_=c(1,dat$s[-length(dat$s)])
  c1t_=c(1,dat$c1t[-length(dat$c1t)])
  c2t_=c(1,dat$c2t[-length(dat$c2t)])
  dt=c(diff(dat$t,1),0)
  ds=dat$s-s_
  w=((n1+n2)*c1t_*c2t_)/((n1*c1t_)+(n2*c2t_))
  dat<-cbind(dat,c1t_,c2t_,dt,w,s_,ds)
  delmax<-dat[1:(length(dat$t)-1),]

  WKM<-with(delmax,sqrt(n1*n2/(n1+n2))*sum((w*(s2t-s1t)*dt)))
  mult<-delmax$w*delmax$s*delmax$dt
  inter<-(cumsum(mult[length(mult):1]))[length(mult):1]
  var22=-sum(with(delmax,(inter^2)*ds/(s*s_*w)))
  Statistic<-WKM/sqrt(var22)
  Pvalue<-2*(1-pnorm(abs(Statistic)))

  result<-data.frame(method=c("Partial Weighted Kaplan-Meier","Partial log-rank"),
                     t0=rep(t0,2),
                     statistic=c(round(Statistic,5),round(logrank.stat,5)),
                     pvalue=c(round(Pvalue,5),round(logrank.P,5)))
  result
}

Short.test<-function(time,status,group,t0){

  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))!=2) stop("There must be two groups")
  if (all(group%in%c(0,1))==FALSE)  stop("Group must be 0 or 1")
  if(t0<=min(time[which(status==1)])) stop("t0 is too small,please select a time point more than the minimum non-censored timepoint")
  if(t0>=max(time[which(status==1)])) stop("t0 is too large,please select a time point less than the maximum of time")

  Short(time,status,group,t0)
}


######################################################################
######         Test homogeneity after a fixed time point        ######
######################################################################
Long<-function(time,status,group,t0){
  group<-group+1
  fit<-survfit(Surv(time,status==1)~group)
  t1<-fit$time[1:fit$strata[[1]]]
  d1<-fit$n.event[1:fit$strata[[1]]]
  Y1<-fit$n.risk[1:fit$strata[[1]]]
  s1t<-fit$surv[1:fit$strata[[1]]]
  n1<-fit$strata[[1]]
  l1<-length(t1[t1<t0])
  H1t<-cumsum(d1/Y1)
  VarH1t<-cumsum(d1/(Y1^2))
  ds1<-data.frame(t1,d1,Y1,s1t,H1t,VarH1t)
  t2<-fit$time[(fit$strata[[1]]+1):(fit$strata[[1]]+fit$strata[[2]])]
  d2<-fit$n.event[(fit$strata[[1]]+1):(fit$strata[[1]]+fit$strata[[2]])]
  Y2<-fit$n.risk[(fit$strata[[1]]+1):(fit$strata[[1]]+fit$strata[[2]])]
  s2t<-fit$surv[(fit$strata[[1]]+1):(fit$strata[[1]]+fit$strata[[2]])]
  n2<-fit$strata[[2]]
  l2<-length(t2[t2<t0])
  H2t<-cumsum(d2/Y2)
  VarH2t<-cumsum(d2/(Y2^2))
  ds2<-data.frame(t2,d2,Y2,s2t,H2t,VarH2t)
  fit1<-survfit(Surv(time,status)~1)
  t<-fit1$time
  n<-n1+n2
  l<-length(t[t<t0])
  d<-fit1$n.event
  Y<-fit1$n.risk
  st2<-fit1$surv^2
  su<-cumsum(d/(Y*(Y-d)))
  VS<-st2*su

  if(t0<max(min(t1),min(t2)))
  {print("t0 is too small,please select another time greater than both minimum time")}
  else
  {H1t0<-ds1$H1t[sum((ds1$t1<=t0)==TRUE)]
  VarH1t0<-ds1$VarH1t[sum((ds1$t1<=t0)==TRUE)]
  H2t0<-ds2$H2t[sum((ds2$t2<=t0)==TRUE)]
  VarH2t0<-ds2$VarH2t[sum((ds2$t2<=t0)==TRUE)]
  XNAt0=H1t0-H2t0
  VarNAt0=VarH1t0+VarH2t0
  ZNAt0<-XNAt0/sqrt(VarNAt0)
  df1<-data.frame(t1,d1,Y1)
  df2<-data.frame(t2,d2,Y2)
  df<-merge(df1,df2,all=TRUE,by.x="t1",by.y="t2")
  df$d1[is.na(df$d1)]<-0
  df$d2[is.na(df$d2)]<-0
  new<-df[which((df$d1!=0)|(df$d2!=0)),]
  tn1<-new[is.na(new$Y1),]$t1
  tn2<-new[is.na(new$Y2),]$t1
  num1<-0
  num2<-0
  for(i in 1:length(tn1))
  {num1[i]<-sum((time[group==1])>=tn1[i])}
  for(j in 1:length(tn2))
  {num2[j]<-sum((time[group==2])>=tn2[j])}
  new[which(is.na(new$Y1)),3]<-num1
  new[which(is.na(new$Y2)),5]<-num2
  if(t0>=max(new$t1))
  {print("prespecified t0 is larger than the max event time in the data")}
  else
  {
    fnset<-subset(new,t1>t0)
    Y<-with(fnset,Y1+Y2)
    d<-with(fnset,d1+d2)
    XLRt0<-with(fnset,sum((Y2*d1-Y1*d2)/Y))
    V<-with(fnset,d*(Y1*Y2/(Y^2))*((Y-d)/(Y-1)))
    V[is.nan(V)]<-0
    VarLRt0<-with(fnset,sum(V))
    ZLRt0<-XLRt0/sqrt(VarLRt0)
  }
  }
  Z1<-ZLRt0
  P1<-2*(1-pnorm(abs(Z1)))
  Z2<-(ZNAt0+ZLRt0)/sqrt(2)
  P2<-2*(1-pnorm(abs(Z2)))
  s1t<-s1t[l1]
  s2t<-s2t[l2]
  VS<-VS[l]
  Z3<-((n1*n2*(s1t-s2t)/n)+XLRt0)/sqrt(n1*n2*VS+VarLRt0)
  P3<-2*(1-pnorm(abs(Z3)))
  Z4<-ZNAt0^2+ZLRt0^2
  P4<-1-pchisq(Z4,df=1)

  result<-data.frame(method=c("Partial log-rank","Zols","Zspp","Qua"),
                     t0=rep(t0,4),
                     statistic=c(round(Z1,5),round(Z2,5),round(Z3,5),round(Z4,5)),
                     pvalue=c(round(P1,5),round(P2,5),round(P3,5),round(P4,5)))
  result
}

Long.test<-function(time,status,group,t0){

  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))!=2) stop("There must be two groups")
  if (all(group%in%c(0,1))==FALSE)  stop("Group must be 0 or 1")
  if(is.na(t0)==TRUE) stop("t0 is null,please select a time point")
  if(t0>=max(time[which(status==1)])) stop("t0 is too large,please select a time point less than the maximum of non-censored time")
  if(t0<=min(time[which(status==1)])) stop("t0 is too small,please select a time point larger than the minimum of non-censored time")

  Long(time,status,group,t0)
}


######################################################################
######               Find for crossing time points              ######
######################################################################
cross<-function(time,status,group){
  dataset<-data.frame(time,status,group)
  all<-survfit(Surv(time,status)~group)

  crosspoint<-c()
  df1<-survfit(Surv(time,status)~1,data=dataset[dataset$group==0,])
  df2<-survfit(Surv(time,status)~1,data=dataset[dataset$group==1,])
  b<-summary(df1)
  c<-summary(df2)
  b$time
  b$surv
  dd<-c()
  tds1<-data.frame(time=b$time,s1t=b$surv)
  tds2<-data.frame(time=c$time,s2t=c$surv)
  tds<-merge(tds1,tds2,all=TRUE)

  for(i in 2:nrow(tds)){
    if(is.na(tds$s1t[i])) tds$s1t[i]<-tds$s1t[i-1]
    if(is.na(tds$s2t[i])) tds$s2t[i]<-tds$s2t[i-1]
  }


  tds$s1t[is.na(tds$s1t)]<-1
  tds$s2t[is.na(tds$s2t)]<-1
  d<-tds$s1t-tds$s2t
  tds<-cbind(tds,d)

  for (i in 2:nrow(tds)){
    dd[i]<-tds$d[i]*tds$d[i-1]
  }
  if (tds$s1t[1]==tds$s2t[1]) dd[1]<-0 else dd[1]<-2
  tds<-cbind(tds,dd)

  for(i in 1:nrow(tds)){
    if (tds$dd[i]<=0) {
      crosspoint1<-tds$time[i]
      crosspoint<-rbind(crosspoint,crosspoint1)
    }
  }
  if (length(crosspoint)==0){
    result=data.frame(crosspoint=NA)} else{
      rownames(crosspoint)<-seq(1:length(crosspoint))
      result<-data.frame(crosspoint)}
  result
}


crosspoint<-function(time,status,group){

  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))!=2) stop("There must be two groups")
  if (all(group%in%c(0,1))==FALSE)  stop("Group must be 0 or 1")

  num<-cross(time,status,group)$crosspoint

  if(sum(is.na(num))) {message("There is no crossing exists.")} else{
    if (length(num)==1) {
      message("\n There is only one crossing exists. \n The exat crossing time point is as follow:\n")
      cross(time,status,group)
    }

    if (length(num)>=2) {
      message("\n There are ", length(num)," crossings exist. \n The exat crossing time points are as follows:\n")
      cross(time,status,group)
    }}
}


########################################################
######            Survival rate plot figure       ######
########################################################

Survival.plot<-function(time,status,group,...){
  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))!=2) stop("There must be two groups")
  if (all(group%in%c(0,1))==FALSE)  stop("Group must be 0 or 1")
  kmfit<-survfit(Surv(time,status==1)~group)
  plot(kmfit,...)
}

########################################################
######     Cumulative hazard rate plot figure     ######
########################################################

Cumhazard.plot<-function(time,status,group,col=c(1,4),lwd=c(1,1),lty=c(1,1),lab.x="",lab.y="",legend=FALSE,local.x=NULL,local.y=NULL,legend.0="",legend.1=""){
  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))!=2) stop("There must be two groups")
  if (all(group%in%c(0,1))==FALSE)  stop("Group must be 0 or 1")

  kmfit<-survfit(Surv(time,status==1)~group)
  s<-summary(kmfit)
  strata<-s$strata
  group0.t<-s$time[which(strata=="group=0")]
  group0.cumhaz<--log(s$surv[which(strata=="group=0")])
  group1.t<-s$time[which(strata=="group=1")]
  group1.cumhaz<--log(s$surv[which(strata=="group=1")])

  plot(c(0,group0.t),c(0,group0.cumhaz),"S",col=col[1],lwd=lwd[1],lty=lty[1],xlab=lab.x,
       ylab=lab.y,cex.lab=1.5,cex.axis=2,
       ylim=c(0,max(group0.cumhaz[is.infinite(group0.cumhaz)==0],group1.cumhaz[is.infinite(group1.cumhaz)==0])),xlim=c(0,max(group0.t,group1.t)))
  lines(c(0,group1.t),c(0,group1.cumhaz),"S",lwd=lwd[2],lty=lty[2],col=col[2])
  if (legend==1) {legend(local.x,local.y,c(legend.0,legend.1),col=col,lwd=lwd,lty=lty,cex=1.5)}
}

########################################################
######      Smooth hazard rate plot figure        ######
########################################################

Hazard.plot<-function(time,status,group,max.0=NULL,max.1=NULL,col=c(1,4),lwd=c(1,1),lty=c(1,1),lab.x="",lab.y="",legend=FALSE,local.x=NULL,local.y=NULL,legend.0="",legend.1=""){
  if (all(status%in%c(0,1))==FALSE) stop("Status must be 0 or 1")
  if (length(unique(group))!=2) stop("There must be two groups")
  if (all(group%in%c(0,1))==FALSE)  stop("Group must be 0 or 1")

  data<-data.frame(time,status,group)
  g0<-data[data$group==0,]
  time0<-g0$time
  status0<-g0$status
  if(is.null(max.0)==TRUE){
    max.0=max(time0)
  }
  fit0<-muhaz(time0,status0,min.time=0,max.time=max.0,bw.method="g",n.est.grid=9,bw.grid=13,kern="b")
  g1<-data[data$group==1,]
  time1<-g1$time
  status1<-g1$status
  if(is.null(max.1)==TRUE){
    max.1=max(time1)
  }
  fit1<-muhaz(time1,status1,min.time=0,max.time=max.1,bw.method="g",n.est.grid=9,bw.grid=13,kern="b")

  haz0<-fit0$haz.est
  haz0.t<-fit0$est.grid
  haz1<-fit1$haz.est
  haz1.t<-fit1$est.grid

  plot(fit0,col=col[1],xlab=lab.x,
       ylab=lab.y,cex.lab=1.5,cex.axis=1.5,lwd=lwd[1],lty=lty[1],
       ylim=c(0,max(haz0,haz1)),xlim=c(0,max(haz0.t,haz1.t)))
  lines(fit1,lwd=lwd[2],lty=lty[2],col=col[2])
  if (legend==1) {legend(local.x,local.y,c(legend.0,legend.1), col=col,lwd=lwd,cex=1.5)}
}

p1<-p2<-1
b1=3
b2=3
n1<-n2<-100
set.seed(1000000*1+1000*n1)
obs.time1<-rexp(n1,1)
obs.status1<-rbinom(n1,size=1,prob=p1)
obs.status1<-ifelse(obs.status1==0,2,1)
cen.time1<-runif(n1,0,b1)
time1<-pmin(obs.time1,cen.time1)
status1<-c(obs.time1<=cen.time1)*obs.status1
group1<-rep(1,n1)
cen.rate1<-length(status1[status1==0])/n1
cause.rate1<-length(status1[status1==1])/n1

set.seed(1000000*2+1000*n2)
obs.time2<-rexp(n2,exp(0.8))
p2<-1-(1-p1)^exp(0.8)
obs.status2<-rbinom(n2,size=1,prob=p2)
obs.status2<-ifelse(obs.status2==0,2,1)
cen.time2<-runif(n2,0,b2)
time2<-pmin(obs.time2,cen.time2)
status2<-c(obs.time2<=cen.time2)*obs.status2
group2<-rep(2,n2)
cen.rate2<-length(status2[status2==0])/n2
cause.rate2<-length(status2[status2==1])/n2



time<-c(time1,time2)
status<-c(status1,status2)
group<-c(group1,group2)-1
PHdata<-data.frame(time,status,group)



b1=6
b2=5

rpweibull<-function (n,A1,A2,t,scale){
  re <- rweibull(n,A1,scale)
  ret1<- re[which(re <= t)]

  ret2<-NULL
  while(length(ret1)<n){
    ind<-n-length(ret1)
    re <- rweibull(ind,A2,scale)
    success <- which(re > t)
    ret2<- re[success]
    ret1<-c(ret1,ret2)
  }
  ret1
}


cen.rate1<-cen.rate2<-c()
cause.rate1<-cause.rate2<-c()
Gray<-R1P<-R2P<-R3P<-R4P<-RP<-c()


set.seed(1000000*1+1000*n1)           #group1:cen
obs.time1<-rpweibull(n1,A1=2,A2=2,t=1,scale=2)
obs.status1<-rbinom(n1,size=1,prob=p1)
obs.status1<-ifelse(obs.status1==0,2,1)
cen.time1<-runif(n1,0,b1)
time1<-pmin(obs.time1,cen.time1)
status1<-c(obs.time1<=cen.time1)*obs.status1
group1<-rep(1,n1)
cen.rate1<-length(status1[status1==0])/n1
cause.rate1<-length(status1[status1==1])/n1

set.seed(1000000*2+1000*n2)           #group2:cen
obs.time2<-rpweibull(n2,A1=0.91,A2=0.91,t=1,scale=2)
obs.status2<-rbinom(n2,size=1,prob=p2)
obs.status2<-ifelse(obs.status2==0,2,1)
cen.time2<-runif(n2,0,b2)
time2<-pmin(obs.time2,cen.time2)
status2<-c(obs.time2<=cen.time2)*obs.status2
group2<-rep(2,n2)
cen.rate2<-length(status2[status2==0])/n2
cause.rate2<-length(status2[status2==1])/n2

time<-c(time1,time2)                                 #final data
status<-c(status1,status2)
group<-c(group1,group2)-1
Crossdata<-data.frame(time,status,group)

#usethis::use_data(Crossdata,PHdata, overwrite = TRUE)

