
#--------------------------------------------------------------------------------------------------------------------------
# Calculating excess hazard and net survival for bias, RMSE and covarage probabilities
#--------------------------------------------------------------------------------------------------------------------------


if(site=="Cervix"){
	
	age.min<-30
	
	modele.true<-as.formula(modele.MI2.cervix)
	cat(paste("\n", "\n", "true modèle","\n",sep=""))
	print(modele.true)
	
	form.true<-as.formula(~I(1/(1+temps)) + I(log(temps+1)) + temps + agex:(I(1/(1+temps)) + I(log(temps+1)) + temps) + agex + 
		I(agex^2) + I(log(agex)) + adiagx + adiagx:temps + (adiagx + 
		adiagx:temps):(I(agex^3) + I((agex^3) * log(agex))))
		
	# ----- > true parameters
	coef.true<-coef.MI2.cervix
	cat(paste("\n", "\n", "true coefs","\n",sep=""))
	print(coef.true)

 }else{
 
	age.min<-45
	
	modele.true<-as.formula(modele.M0.oeso)
	cat(paste("\n", "\n", "true modèle","\n",sep=""))
	print(modele.true)
	
	form.true<-as.formula(~I(1/(1+temps)) + I(log(temps+1)) + temps + agex:(I(1/(1+temps)) + I(log(temps+1)) + temps) + agex + 
		I(agex^-0.5) + I(log(agex)))
		
	# ----- > true parameters
	coef.true<-coef.M0.oeso
	cat(paste("\n", "\n", "true coefs","\n",sep=""))
	print(coef.true)
 
 }



 
sink(paste(Resultats,"/log.calculs.",site,".",n_data,".txt",sep=""))



print(system("hostname", intern=T))
heure.debut=Sys.time()
print(heure.debut)

sessionInfo()


load(paste(Resultats,"/ajust.",site,".",n_data,".RData",sep=""))


######## Predictions ----------------------------------------------

newt<-c(0.5,1,2,3,4,5)
len<-length(newt)

cat("\n","\n","Valeurs de prédictions","\n")


print(newt)
cat("\n","\n","-------------------","\n","\n")

age<-seq(age.min,85,by=5)
print(age)
cat("\n","\n","-------------------","\n","\n")

annee<-seq(1990,2010,by=5)
print(annee)
cat("\n","\n","-------------------","\n","\n")

len.age<-length(age)
len.annee<-length(annee)

 

# Hazard
pred.haz<-function(formule,beta1,VarF,VarB,VarC,Z.tensor,list.tensor,newdata,conf.int=0.95){

	qt.norm<-qnorm(1-(1-conf.int)/2)

	myMat<-design.matrix(formule,data.spec=newdata,Z.add=NULL,Z.tensor=Z.tensor,Z.tint=NULL,list.add=NULL,list.tensor=list.tensor,list.tint=NULL)

	pred<-as.vector(myMat%*%beta1)
	haz<-as.vector(exp(pred))
	
	
	if (!is.null(VarF)){
		se.freq<-sqrt(rowSums((myMat%*%VarF)*myMat))
		haz.freq.inf<-as.vector(exp(pred-qt.norm*se.freq))
		haz.freq.sup<-as.vector(exp(pred+qt.norm*se.freq))
	}else{
		se.freq=NULL
		haz.freq.inf<-NULL
		haz.freq.sup<-NULL
	}
	
	if (!is.null(VarB)){
		se.bayes<-sqrt(rowSums((myMat%*%VarB)*myMat))
		haz.inf<-as.vector(exp(pred-qt.norm*se.bayes))
		haz.sup<-as.vector(exp(pred+qt.norm*se.bayes))
	}else{
		se.bayes=NULL
		haz.inf<-NULL
		haz.sup<-NULL
	}
	
	if (!is.null(VarC)){

		se.cor<-sqrt(rowSums((myMat%*%VarC)*myMat))
		haz.inf.cor<-as.vector(exp(pred-qt.norm*se.cor))
		haz.sup.cor<-as.vector(exp(pred+qt.norm*se.cor))

	}else{
		se.cor=NULL
		haz.inf.cor<-NULL
		haz.sup.cor<-NULL

	}

	list(haz=haz,haz.freq.inf=haz.freq.inf,haz.freq.sup=haz.freq.sup,
	haz.inf=haz.inf,haz.sup=haz.sup,haz.inf.cor=haz.inf.cor,haz.sup.cor=haz.sup.cor,
	se.bayes=se.bayes,se.freq=se.freq,se.cor=se.cor)
}

# Survival
pred.surv <-function(formule,beta1,VarF,VarB,VarC,Z.tensor,list.tensor,newdata,n.legendre=50,conf.int=0.95){

	qt.norm<-qnorm(1-(1-conf.int)/2)

	t1<-newdata[,"temps"]

	t0<-rep(0,length(t1))

	tm<-(t1-t0)/2

	leg<-statmod::gauss.quad(n=n.legendre,kind="legendre")

	# Design matrices for Gauss-Legendre quadrature
	X.func<-function(t1,data,formule,Z.tensor,list.tensor){

		data.t<-data
		data.t[,"temps"]<-t1
		design.matrix(formule,data.spec=data.t,Z.add=NULL,Z.tensor=Z.tensor,Z.tint=NULL,list.add=NULL,list.tensor=list.tensor,list.tint=NULL)

	}

	X.GL<-lapply(1:n.legendre, function(i) X.func(tm*leg$nodes[i]+(t0+t1)/2,newdata,formule,Z.tensor,list.tensor))

	cumul.haz <- lapply(1:n.legendre, function(i) as.vector(exp((X.GL[[i]]%*%beta1)))*leg$weights[i])

	cumul.haz<-tm*Reduce("+",cumul.haz)

	surv=exp(-cumul.haz)
	
	# if any cumul hazard is zero, we put it at a very low positive value
	cumul.haz[cumul.haz==0]<-1e-16

	deriv.cumul.haz <- lapply(1:n.legendre, function(i) X.GL[[i]]*as.vector(exp((X.GL[[i]]%*%beta1)))*leg$weights[i])

	deriv.cumul.haz<-tm*Reduce("+",deriv.cumul.haz)

	#---------------------------
	log.cumul.haz<-log(cumul.haz)

	deriv.log.cumul.haz<-deriv.cumul.haz/cumul.haz

	#---------------------------
	# Delta method
		
	if (!is.null(VarB)){
	
		se.bayes<-sqrt(rowSums((deriv.log.cumul.haz%*%VarB)*deriv.log.cumul.haz))
		surv.inf=exp(-exp(log.cumul.haz-qt.norm*se.bayes))
		surv.sup=exp(-exp(log.cumul.haz+qt.norm*se.bayes))
		
	}else{
		se.bayes=NULL
		surv.inf=NULL
        surv.sup=NULL

	}	
		
	if (!is.null(VarF)){	

		se.freq<-sqrt(rowSums((deriv.log.cumul.haz%*%VarF)*deriv.log.cumul.haz))
		surv.freq.inf=exp(-exp(log.cumul.haz-qt.norm*se.freq))
		surv.freq.sup=exp(-exp(log.cumul.haz+qt.norm*se.freq))
	}else{	
		se.freq=NULL
		surv.freq.inf=NULL
		surv.freq.sup=NULL
	}
	
		
	if (!is.null(VarC)){

		se.cor<-sqrt(rowSums((deriv.log.cumul.haz%*%VarC)*deriv.log.cumul.haz))
		surv.inf.cor=exp(-exp(log.cumul.haz-qt.norm*se.cor))
		surv.sup.cor=exp(-exp(log.cumul.haz+qt.norm*se.cor))

	}else{
		se.cor=NULL
		surv.inf.cor=NULL
        surv.sup.cor=NULL

	}
	#---------------------------
	list(surv=surv,surv.freq.inf=surv.freq.inf,surv.freq.sup=surv.freq.sup,surv.inf=surv.inf,surv.sup=surv.sup,
	surv.inf.cor=surv.inf.cor,surv.sup.cor=surv.sup.cor,se.bayes=se.bayes,se.freq=se.freq,se.cor=se.cor)

}


#---------- functions to predict survival with the true model
Rescale.in <- function(gl,a,b){
  gl$nodes <- gl$nodes*(b-a)/2+(a+b)/2
  gl$weights <- gl$weights*(b-a)/2
  return(gl)
	}
survieFunc.in<-function(x,beta, modele,GL,a=0){
  gg=Rescale.in(GL,a,x["times"])
  myMat=model.matrix(as.formula(modele),data.frame(intnum=gg$nodes,unsurt=1/(gg$nodes+1),logt=log(gg$nodes+1),agex=as.numeric(x["agex"]),adiagx=as.numeric(x["adiagx"])))
  return(exp(-sum(exp(myMat%*%beta)*gg$weights)))
	}
pred.surv.in <- function(  times , agex , adiagx , beta, modele)
	{
	GL=statmod::gauss.quad(n=50,kind="legendre")
	temp <- data.frame( times = rep( times , length(agex * length(adiagx) ))         ,
                      agex  = rep( rep( agex , each = length(times) ) , length(adiagx) )          ,
                      adiagx= rep( adiagx , each = length(times) * length(agex) )                 )                     
	
	temp$surv <- apply( temp , 1 , FUN =survieFunc.in, beta=beta, modele=modele, GL=GL  )
	return(temp$surv)
	}

#----------------------------------------------------------------------------


# indices of files which present convergence issues

true.conv<-fail.true
LCV.conv<-fail.LCV
LAML.conv<-fail.LAML

res.simu<-data.frame(
time=rep(newt,len.age*len.annee*NbFichier),
age=rep(age,each=len),
annee=rep(annee,each=len.age*len),
fichier=rep(1:NbFichier,each=len*len.age*len.annee))

# convergence indicator
res.simu$conv.true=1
res.simu$conv.LCV=1
res.simu$conv.LAML=1
if(length(true.conv)!=0) res.simu[res.simu$fichier%in%true.conv,]$conv.true=0
if(length(LCV.conv)!=0) res.simu[res.simu$fichier%in%LCV.conv,]$conv.LCV=0
if(length(LAML.conv)!=0) res.simu[res.simu$fichier%in%LAML.conv,]$conv.LAML=0

res.simu$haz.true<-0
res.simu$surv.true<-0
res.simu$haz.fit.true<-0
res.simu$haz.LCV<-0
res.simu$haz.LAML<-0
#----------------------------
res.simu$haz.fit.true.se<-0
res.simu$haz.LCV.se<-0
res.simu$haz.LAML.se<-0
#----------------------------
res.simu$haz.fit.true.se.freq<-0
res.simu$haz.LCV.se.freq<-0
res.simu$haz.LAML.se.freq<-0
#------------------------------
res.simu$haz.LAML.se.cor<-0
#------------------------------
res.simu$surv.fit.true<-0
res.simu$surv.LCV<-0
res.simu$surv.LAML<-0
#------------------------------
res.simu$surv.fit.true.se<-0
res.simu$surv.LCV.se<-0
res.simu$surv.LAML.se<-0
#------------------------------
res.simu$surv.fit.true.se.freq<-0
res.simu$surv.LCV.se.freq<-0
res.simu$surv.LAML.se.freq<-0
#------------------------------
res.simu$surv.LAML.se.cor<-0


long<-len*len.age*len.annee



for (age_k in age){

	for (annee_k in annee){

		# true model
		newdata=data.frame(temps=newt,agec=trunc(age_k)-70,adiagc=trunc(annee_k)-2000,unsurt=1/(1+newt),logt=log(newt+1),intnum=newt,agex=trunc(age_k)/10,adiagx=trunc(annee_k)-1989)
				
		X.true<-model.matrix(modele.true,newdata)
					
		haz.true<-as.vector(exp(X.true%*%coef.true))

		surv.true<-pred.surv.in(newt,trunc(age_k)/10,trunc(annee_k)-1989,coef.true,modele.true)
		
		res.simu[(res.simu$age==age_k)&(res.simu$annee==annee_k),]$haz.true<-rep(haz.true,NbFichier)
		res.simu[(res.simu$age==age_k)&(res.simu$annee==annee_k),]$surv.true<-rep(surv.true,NbFichier)
		
		
		# predictions with models LAML and lCV

		for (j in 1:NbFichier){

			res<-res.tot[[j]]
			
			res.true<-res$res.true
			res.LCV<-res$res.LCV
			res.LAML<-res$res.LAML
			
			k.t<-res$k.t
			k.agec<-res$k.agec
			k.adiagc<-res$k.adiagc
			
			formule<-as.formula(~tensor(temps,agec,adiagc,knots=list(k.t, k.agec,k.adiagc)))
			
			
			if(length(res.true)!=1){
				predit.haz.true<-pred.haz(form.true,res.true$beta.true,res.true$VarF.true,res.true$VarB.true,NULL,NULL,NULL,newdata)
				predit.surv.true<-pred.surv(form.true,res.true$beta.true,res.true$VarF.true,res.true$VarB.true,NULL,NULL,NULL,newdata)
			}else{
				predit.haz.true<-NULL
				predit.surv.true<-NULL
			}
				
			if(length(res.LCV)!=1){
				predit.haz.LCV<-pred.haz(formule,res.LCV$beta.LCV,res.LCV$VarF.LCV,res.LCV$VarB.LCV,NULL,res.LCV$Z.tensor,res.LCV$list.tensor,newdata)
				predit.surv.LCV<-pred.surv(formule,res.LCV$beta.LCV,res.LCV$VarF.LCV,res.LCV$VarB.LCV,NULL,res.LCV$Z.tensor,res.LCV$list.tensor,newdata)
			}else{
				predit.haz.LCV<-NULL
				predit.surv.LCV<-NULL
			}
				
			if(length(res.LAML)!=1){
				predit.haz.LAML<-pred.haz(formule,res.LAML$beta.LAML,res.LAML$VarF.LAML,res.LAML$VarB.LAML,res.LAML$VarC,res.LAML$Z.tensor,res.LAML$list.tensor,newdata)
				predit.surv.LAML<-pred.surv(formule,res.LAML$beta.LAML,res.LAML$VarF.LAML,res.LAML$VarB.LAML,res.LAML$VarC,res.LAML$Z.tensor,res.LAML$list.tensor,newdata)
			}else{
				predit.haz.LAML<-NULL
				predit.surv.LAML<-NULL
			}	
				
				
			condition<-(res.simu$age==age_k)&(res.simu$annee==annee_k)&(res.simu$fichier==j)
			# hazard
			res.simu[condition,]$haz.fit.true<-predit.haz.true$haz
			res.simu[condition,]$haz.LCV<-predit.haz.LCV$haz
			res.simu[condition,]$haz.LAML<-predit.haz.LAML$haz

			# Frequentist standard errors for hazard
			res.simu[condition,]$haz.fit.true.se.freq<-predit.haz.true$se.freq
			res.simu[condition,]$haz.LCV.se.freq<-predit.haz.LCV$se.freq
			res.simu[condition,]$haz.LAML.se.freq<-predit.haz.LAML$se.freq
			
			# Bayesian standard errors for hazard
			res.simu[condition,]$haz.fit.true.se<-predit.haz.true$se.bayes
			res.simu[condition,]$haz.LCV.se<-predit.haz.LCV$se.bayes
			res.simu[condition,]$haz.LAML.se<-predit.haz.LAML$se.bayes
			
			# Corrected bayesian standard errors for hazard
			res.simu[condition,]$haz.LAML.se.cor<-predit.haz.LAML$se.cor

			# survival
			res.simu[condition,]$surv.fit.true<-predit.surv.true$surv
			res.simu[condition,]$surv.LCV<-predit.surv.LCV$surv
			res.simu[condition,]$surv.LAML<-predit.surv.LAML$surv
			
			# Frequentist Confidence intervals for survival
			res.simu[condition,]$surv.fit.true.se.freq<-predit.surv.true$se.freq
			res.simu[condition,]$surv.LCV.se.freq<-predit.surv.LCV$se.freq
			res.simu[condition,]$surv.LAML.se.freq<-predit.surv.LAML$se.freq
			
			# Bayesian Confidence intervals for survival
			res.simu[condition,]$surv.fit.true.se<-predit.surv.true$se.bayes
			res.simu[condition,]$surv.LCV.se<-predit.surv.LCV$se.bayes
			res.simu[condition,]$surv.LAML.se<-predit.surv.LAML$se.bayes
			
			# Corrected bayesian standard errors for survival
			res.simu[condition,]$surv.LAML.se.cor<-predit.surv.LAML$se.cor
			
		}

		
		cat("\n","\n","age ",age_k," et annee",annee_k," ok","\n","\n")
		
		
	}
}

	res.simu$elapsed.true<-0
	res.simu$elapsed.LCV<-0
	res.simu$elapsed.LAML<-0
	res.simu$ll.true<-0
	res.simu$ll.unpen.LCV<-0
	res.simu$ll.unpen.LAML<-0
	res.simu$ll.LCV<-0
	res.simu$ll.LAML<-0
	res.simu$max.grad.beta.true<-0
	res.simu$max.grad.beta.LCV<-0
	res.simu$max.grad.beta.LAML<-0
	res.simu$max.grad.lambda.LCV<-0
	res.simu$max.grad.lambda.LAML<-0
	res.simu$aic.true<-0
	res.simu$aic.LCV<-0
	res.simu$aic.LAML<-0
	res.simu$aicc<-0

for (j in 1:NbFichier){

	res<-res.tot[[j]]
		
	res.true<-res$res.true
	res.LCV<-res$res.LCV
	res.LAML<-res$res.LAML
	
	if(length(res.true)!=1){
	
		res.simu$elapsed.true[(1+(j-1)*long):(j*long)]<-rep(res.true$elapsed.true,long)
		res.simu$ll.true[(1+(j-1)*long):(j*long)]<-rep(res.true$ll.true,long)
		
		res.simu$max.grad.beta.true[(1+(j-1)*long):(j*long)]<-rep(max(abs(res.true$grad.beta.true)),long)
		
		res.simu$aic.true[(1+(j-1)*long):(j*long)]<-rep(res.true$aic.true,long)
		
	}else{
	
		res.simu$elapsed.true[(1+(j-1)*long):(j*long)]<-0
		res.simu$ll.true[(1+(j-1)*long):(j*long)]<-0
		
		res.simu$max.grad.beta.true[(1+(j-1)*long):(j*long)]<-0
		
		res.simu$aic.true[(1+(j-1)*long):(j*long)]<-0
		
	}
	
	if(length(res.LCV)!=1){
	
		res.simu$elapsed.LCV[(1+(j-1)*long):(j*long)]<-rep(res.LCV$elapsed.LCV,long)
		res.simu$ll.unpen.LCV[(1+(j-1)*long):(j*long)]<-rep(res.LCV$ll.unpen.LCV,long)
		res.simu$ll.LCV[(1+(j-1)*long):(j*long)]<-rep(res.LCV$ll.LCV,long)
		
		res.simu$max.grad.beta.LCV[(1+(j-1)*long):(j*long)]<-rep(max(abs(res.LCV$grad.beta.LCV)),long)
		res.simu$max.grad.lambda.LCV[(1+(j-1)*long):(j*long)]<-rep(max(abs(res.LCV$grad.lambda.LCV)),long)
		
		res.simu$aic.LCV[(1+(j-1)*long):(j*long)]<-rep(res.LCV$aic.LCV,long)

	}else{
	
		res.simu$elapsed.LCV[(1+(j-1)*long):(j*long)]<-0
		res.simu$ll.unpen.LCV[(1+(j-1)*long):(j*long)]<-0
		res.simu$ll.LCV[(1+(j-1)*long):(j*long)]<-0
	
		res.simu$max.grad.beta.LCV[(1+(j-1)*long):(j*long)]<-0
		res.simu$max.grad.lambda.LCV[(1+(j-1)*long):(j*long)]<-0
		
		res.simu$aic.LCV[(1+(j-1)*long):(j*long)]<-0

	}
	
	if(length(res.LAML)!=1){
	
		res.simu$elapsed.LAML[(1+(j-1)*long):(j*long)]<-rep(res.LAML$elapsed.LAML,long)
		res.simu$ll.unpen.LAML[(1+(j-1)*long):(j*long)]<-rep(res.LAML$ll.unpen.LAML,long)
		res.simu$ll.LAML[(1+(j-1)*long):(j*long)]<-rep(res.LAML$ll.LAML,long)
		
		res.simu$max.grad.beta.LAML[(1+(j-1)*long):(j*long)]<-rep(max(abs(res.LAML$grad.beta.LAML)),long)
		res.simu$max.grad.lambda.LAML[(1+(j-1)*long):(j*long)]<-rep(max(abs(res.LAML$grad.lambda.LAML)),long)
		
		res.simu$aic.LAML[(1+(j-1)*long):(j*long)]<-rep(res.LAML$aic.LAML,long)
		res.simu$aicc[(1+(j-1)*long):(j*long)]<-rep(res.LAML$aicc,long)

	}else{
	
		res.simu$elapsed.LAML[(1+(j-1)*long):(j*long)]<-0
		res.simu$ll.unpen.LAML[(1+(j-1)*long):(j*long)]<-0
		res.simu$ll.LAML[(1+(j-1)*long):(j*long)]<-0
	
		res.simu$max.grad.beta.LAML[(1+(j-1)*long):(j*long)]<-0
		res.simu$max.grad.lambda.LAML[(1+(j-1)*long):(j*long)]<-0
		
		res.simu$aic.LAML[(1+(j-1)*long):(j*long)]<-0
		res.simu$aicc[(1+(j-1)*long):(j*long)]<-0

	}
		
}


save(res.simu,file = paste(Resultats,"/res.taux.surv.",site,".",n_data,".RData",sep="") ,compress="xz")

rm(res.simu)



cat("--------------------------------------------","\n")
cat("--------------------------------------------","\n")

heure.fin=Sys.time()
cat("\n", "\n")
print(heure.fin)

cat("\n", "\n", "programme's duration","\n")
print(difftime(heure.fin, heure.debut,units = c("hours")))



sink()

########################################






