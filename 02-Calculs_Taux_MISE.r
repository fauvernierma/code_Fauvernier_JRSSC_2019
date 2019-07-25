


#--------------------------------------------------------------------------------------------------------------------------
# Calculating excess hazard for RMISE (using GL quadrature with 20 points) 
#--------------------------------------------------------------------------------------------------------------------------


sink(paste(Resultats,"/log.",site,".",n_data,".txt",sep=""),append=TRUE)
load(paste(Resultats,"/ajust.",site,".",n_data,".RData",sep=""))


######## Predictions ----------------------------------------------

# Hazard
pred.haz<-function(formule,beta1,Z.tensor,list.tensor,newdata,conf.int=0.95){

	qt.norm<-qnorm(1-(1-conf.int)/2)

	myMat<-design.matrix(formule,data.spec=newdata,Z.add=NULL,Z.tensor=Z.tensor,Z.tint=NULL,list.add=NULL,list.tensor=list.tensor,list.tint=NULL)

	pred<-as.vector(myMat%*%beta1)
	haz<-as.vector(exp(pred))
	
	haz
}


#---------- to adapt points and weights for GL quadrature on the interval [0;5]
Rescale.in <- function(gl,a,b){
  gl$nodes <- gl$nodes*(b-a)/2+(a+b)/2
  gl$weights <- gl$weights*(b-a)/2
  return(gl)
	}

#----------------------------------------------------------------------------

GL<-statmod::gauss.quad(20)

GL5<-Rescale.in(GL,0,5)

newt<-GL5$nodes
len<-length(newt)

age<-seq(age.min,90,by=5)
annee<-seq(1990,2010,by=5)

len.age<-length(age)
len.annee<-length(annee)


# indices of files which present convergence issues
true.conv<-fail.true
if(length(true.conv)!=0) {
	cat("indices : true models that did not converge","\n")
	print(true.conv)
}

LCV.conv<-fail.LCV
if(length(LCV.conv)!=0) {
	cat("indices : LCV models that did not converge","\n")
	print(LCV.conv)
}

LAML.conv<-fail.LAML
if(length(LAML.conv)!=0) {
	cat("indices : LAML models that did not converge","\n")
	print(LAML.conv)
}
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
		
long<-len*len.age*len.annee


for (age_k in age){

	for (annee_k in annee){

		# true model
		newdata=data.frame(temps=newt,agec=trunc(age_k)-70,adiagc=trunc(annee_k)-2000,unsurt=1/(1+newt),logt=log(newt+1),intnum=newt,agex=trunc(age_k)/10,adiagx=trunc(annee_k)-1989)
				
		X.true<-model.matrix(modele.true,newdata)
					
		haz.true<-as.vector(exp(X.true%*%coef.true))

		res.simu[(res.simu$age==age_k)&(res.simu$annee==annee_k),]$haz.true<-rep(haz.true,NbFichier)
		
		# predictions with LAML and LCV models

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
				predit.haz.true<-pred.haz(form.true,res.true$beta.true,NULL,NULL,newdata)
			}else{
				predit.haz.true<-NULL
			}
				
			if(length(res.LCV)!=1){
				predit.haz.LCV<-pred.haz(formule,res.LCV$beta.LCV,res.LCV$Z.tensor,res.LCV$list.tensor,newdata)
			}else{
				predit.haz.LCV<-NULL
			}
				
			if(length(res.LAML)!=1){
				predit.haz.LAML<-pred.haz(formule,res.LAML$beta.LAML,res.LAML$Z.tensor,res.LAML$list.tensor,newdata)
			}else{
				predit.haz.LAML<-NULL
			}	
				
				
			condition<-(res.simu$age==age_k)&(res.simu$annee==annee_k)&(res.simu$fichier==j)
			# hazard
			res.simu[condition,]$haz.fit.true<-predit.haz.true
			res.simu[condition,]$haz.LCV<-predit.haz.LCV
			res.simu[condition,]$haz.LAML<-predit.haz.LAML

		}

	}
}

	
save(res.simu,file = paste(Resultats,"/res.MISE.",site,".",n_data,".RData",sep="") ,compress="xz")

rm(res.simu)
sink()

########################################






