
#--------------------------------------------------------------------------------------------------------------------------
# Fiting all models 
#--------------------------------------------------------------------------------------------------------------------------

sink(paste(Resultats,"/log.",site,".",n_data,".txt",sep=""),append=TRUE)
res.tot<-vector("list",NbFichier)

fail.true <- NULL
fail.LCV <- NULL
fail.LAML <- NULL

for (i in 1:NbFichier){

	tempdata=ListDataSimulation[[i]]
	#----------------------------------------------------------------------------------------

	# ----- > adding necessary information
	tempdata$etat5=ifelse(tempdata$DureeDeVie<5,tempdata$etat,0)
	tempdata$DureeDeVie5=ifelse(tempdata$DureeDeVie<5,tempdata$DureeDeVie,5)

	jeudata<-tempdata
	jeudata$anneec  <- jeudata$annee.entier - 1960

	jeudata$age.exit <-  jeudata$age.entier + jeudata$DureeDeVie
	jeudata$age.exit <-  ifelse(jeudata$age.exit >= 99, 99, jeudata$age.exit)

	jeudata$age.exit.r <- floor(jeudata$age.exit)

	jeudata$annee.exit <-  jeudata$annee.entier + jeudata$DureeDeVie
	
	jeudata$annee.exit.r <- floor(jeudata$annee.exit)

	jeudata <- merge(jeudata,tauxatt,by.x=c("annee.exit.r","age.exit.r","sexe"),by.y=c("ANNEE","AGEX","SEXE"))

	jeudata$agec=trunc(jeudata$age)-70
	jeudata$adiagc=trunc(jeudata$annee)-2000


	jeudata$temps<-jeudata$DureeDeVie5
	jeudata$Fail<-jeudata$etat5

	cat("\n")
	print(summary(jeudata$temps))
	cat("\n")	

	jeudata$unsurt=1/(1+jeudata$temps)
	jeudata$logt=log(jeudata$temps+1)
	jeudata$intnum=jeudata$temps
	jeudata$agex=trunc(jeudata$age)/10
	jeudata$adiagx=trunc(jeudata$annee)-1989
		
	######## knots ----------------------------------------------
	k.t=quantile(tempdata[tempdata$etat5==1,c("DureeDeVie5")], p=c(0, 0.20, 0.40, 0.60, 0.80, 1))
	cat("\n", "\n", ".......... knots for time.......","\n")
	print(k.t)


	k.agec=-70+quantile(tempdata[tempdata$etat5==1,c("age.entier")], p=c(0, 0.25, 0.50, 0.75, 1))
	cat("\n", "\n", ".......... knots for age.......","\n")
	print(k.agec)

	k.adiagc=-2000+quantile(tempdata[tempdata$etat5==1,c("annee.entier")], p=c(0, 0.33, 0.66, 1))
	cat("\n", "\n", ".......... knots for year.......","\n")
	print(k.adiagc)


	# fitting the true model
	time.true<-system.time(mod.true<-try(survPen(form.true,jeudata,t1=temps,event=Fail,expected=MUA,type="net",criterion="LAML",tol.beta=tol.beta1,tol.rho=tol.sp1,max.it.beta=max.it1)))
	
	if (class(mod.true)[1]!="survPen") {
		mod.true$converged=FALSE
		cat(paste("\n", "\n", "mod.true ",i," error","\n",sep=""))
	}else{
		if(mod.true$converged){
		cat(paste("\n", "\n", "mod.true ",i," ok","\n",sep=""))
		
		
		}else{
		cat(paste("\n", "\n", "mod.true ",i," does not converge","\n",sep=""))
		}

	}

	
	
	# tensor formula
	form.simu<- ~tensor(temps,agec,adiagc,knots=list(k.t, k.agec,k.adiagc))

	
	# fitting tensor with LCV
	time.LCV<-system.time(mod.LCV<-try(survPen(form.simu,jeudata,t1=temps,event=Fail,expected=MUA,type="net",criterion="LCV",tol.beta=tol.beta1,tol.rho=tol.sp1,max.it.beta=max.it1,max.it.rho=max.it2)))
	
	if (class(mod.LCV)[1]!="survPen") {
		mod.LCV$converged=FALSE
		cat(paste("\n", "\n", "mod.LCV ",i," error","\n",sep=""))
	}else{
		if(mod.LCV$converged){
		cat(paste("\n", "\n", "mod.LCV ",i," ok","\n",sep=""))
		
		}else{
		cat(paste("\n", "\n", "mod.LCV ",i," does not converge","\n",sep=""))
		}
	}
	
	# fitting tensor with LAML
	time.LAML<-system.time(mod.LAML<-try(survPen(form.simu,jeudata,t1=temps,event=Fail,expected=MUA,type="net",criterion="LAML",tol.beta=tol.beta1,tol.rho=tol.sp1,max.it.beta=max.it1,max.it.rho=max.it2)))
	
	if (class(mod.LAML)[1]!="survPen") {
		mod.LAML$converged=FALSE
		cat(paste("\n", "\n", "mod.LAML ",i," error","\n",sep=""))
	}else{
		if(mod.LAML$converged){
		cat(paste("\n", "\n", "mod.LAML ",i," ok","\n",sep=""))
		
		
		}else{
		cat(paste("\n", "\n", "mod.LAML ",i," does not converge","\n",sep=""))
		}
	}
	
	
	if (mod.true$converged){
		res.true<-list(elapsed.true=round(time.true[3],0),beta.true=mod.true$beta.hat,
		grad.beta.true=mod.true$grad.beta,VarF.true=mod.true$VarF,VarB.true=mod.true$VarB,
		ll.true=mod.true$ll,aic.true=mod.true$aic,edf=mod.true$edf)
	}else{
	
		res.true<-"fail"
		fail.true <- c(fail.true,i)
	
	}

	
	if (mod.LCV$converged){
		res.LCV<-list(elapsed.LCV=round(time.LCV[3],0),list.tensor=mod.LCV$list.tensor,
		Z.tensor=mod.LCV$Z.tensor,beta.LCV=mod.LCV$beta.hat,lambda.LCV=mod.LCV$lambda,
		grad.beta.LCV=mod.LCV$grad.beta,grad.lambda.LCV=mod.LCV$grad.rho,
		VarF.LCV=mod.LCV$VarF,VarB.LCV=mod.LCV$VarB,
		LCV=mod.LCV$LCV,ll.LCV=mod.LCV$ll,ll.unpen.LCV=mod.LCV$ll.unpen,
		aic.LCV=mod.LCV$aic,edf=mod.LCV$edf)
	}else{
	
		res.LCV<-"fail"
		fail.LCV <- c(fail.LCV,i)
	
	}
	
	if (mod.LAML$converged){
		res.LAML<-list(elapsed.LAML=round(time.LAML[3],0),list.tensor=mod.LAML$list.tensor,
		Z.tensor=mod.LAML$Z.tensor,beta.LAML=mod.LAML$beta.hat,lambda.LAML=mod.LAML$lambda,
		grad.beta.LAML=mod.LAML$grad.beta,grad.lambda.LAML=mod.LAML$grad.rho,
		VarF.LAML=mod.LAML$VarF,VarB.LAML=mod.LAML$VarB,
		VarC=mod.LAML$VarC,LAML=mod.LAML$LAML,ll.LAML=mod.LAML$ll,
		ll.unpen.LAML=mod.LAML$ll.unpen,aic.LAML=mod.LAML$aic,aicc=mod.LAML$aicc,edf=mod.LAML$edf)
	}else{
	
		res.LAML<-"fail"
		fail.LAML <- c(fail.LAML,i)
	
	}
	
	res.tot[[i]]<-list(res.true=res.true,res.LCV=res.LCV,res.LAML=res.LAML,
	k.t=k.t, k.agec=k.agec,k.adiagc=k.adiagc)
	
	
}	


save(fail.true,fail.LCV,fail.LAML,res.tot,file = paste(Resultats,"/ajust.",site,".",n_data,".RData",sep="") ,compress="xz")

rm(res.tot)
sink()

#############################################################################


