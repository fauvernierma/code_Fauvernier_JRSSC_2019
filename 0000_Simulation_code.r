

######################################################################################################
######################################################################################################
######################################################################################################

# This file contains code used in the simulation study 

#--------------------------------------------------------------------------------------------------------------------------
# Prerequisites
#--------------------------------------------------------------------------------------------------------------------------

# 1)	Need R>=3.5.0
# change this path (finish by /)
Dir <- "path/" 

# for the results
Resultats<-paste0(Dir,"Resultats")
setwd(Resultats)  
 
 

# load survPen v0.0.1.100 used for the article
install.packages(paste0(Dir,"survPen_0.0.1.1000.tar.gz"),repos = NULL, type="source")
library(survPen)
 

# for comparison with rstpm2 (each model fitted with rstpm2 takes about 800 seconds to run)
# library(rstpm2)  # rstpm2 v1.4.5 recommended
# in order to fit rstpm2 models (see Supplementary Material C.2), the following code was used :

#df3 <- c(6,5,4) # for rstpm2, we use the default knots (this is of minor importance when it comes to the simulation results)
	
#	time.rstpm2<-system.time(mod.rstpm2<-try(pstpm2(Surv(temps,Fail)~1, data=jeudata, smooth.formula = ~te(log(temps),agec,adiagc,bs="cr",k=df3),
#		tvc = NULL,
#		bhazard = jeudata$MUA,
#		sp=NULL,
#		criterion="GCV",
#		link.type="PH")))

# to predict using newdata with rstpm2, use the following code
# pred.haz <- predict(mod.rstpm2,newdata,type="haz")
# pred.surv <- predict(mod.rstpm2,newdata,type="surv")


# choose your scenario and sample size
# site is "Oeso" or "Cervix"
# n_data is "2000", "5000" or "10000"

site<-"Oeso"
n_data<-"2000"

# Be careful, if you want to fit all the files 
# this may take between 2 and 4 days (per scenario * sample size) 
if(n_data=="10000"){
	
	NbFichier=200

 }else{
 
	NbFichier=1000
 
 }

# survPen options regarding convergence criteria and maximum number of iterations
tol.beta1=1e-04
tol.sp1=1e-04
max.it1=200
max.it2=30

sink(paste(Resultats,"/log.",site,".",n_data,".txt",sep=""))

print(system("hostname", intern=T))
heure.debut=Sys.time()
print(heure.debut)

sessionInfo()

cat(paste("\n", "\n", "Site","\n",sep=""))
print(site)
cat(paste("\n", "\n", "Nombre d'individus","\n",sep=""))
print(n_data)



#---------------------------------------------------------------------------------------
# LOADING DATA
#---------------------------------------------------------------------------------------


# ----- > expected mortality rates
tauxatt=dget(file=paste0(Dir,"muaDF8933.dat")  )


if(site=="Cervix"){
	
	# age
	load(file =paste0(Dir,"age.cervix.2000.Rdata" ))
	
	age.min<-round(quantile(age.cervix.2000,0.05)/5,0)*5
	cat(paste("\n", "\n", "minimum age","\n",sep=""))
	print(age.min)
	# ----- > loading formula for true model
	load(file =paste0(Dir,"modele.MI2.cervix.Rdata" ) )
	modele.true<-as.formula(modele.MI2.cervix)
	cat(paste("\n", "\n", "true model","\n",sep=""))
	print(modele.true)
	
	form.true<-as.formula(~I(1/(1+temps)) + I(log(temps+1)) + temps + agex:(I(1/(1+temps)) + I(log(temps+1)) + temps) + agex + 
		I(agex^2) + I(log(agex)) + adiagx + adiagx:temps + (adiagx + 
		adiagx:temps):(I(agex^3) + I((agex^3) * log(agex))))
		
	# ----- > loading true parameters
	load(file =paste0(Dir,"coef.MI2.cervix.Rdata" ))
	coef.true<-coef.MI2.cervix
	cat(paste("\n", "\n", "true coefs","\n",sep=""))
	print(coef.true)
	
	if(n_data=="2000"){
		# ----- > loading simulated datasets
		load(paste0(Dir,"ListDataSimulation.MI2.RData"))
		load(paste0(Dir,"DataDesign.MI2.RData"))

		ListDataSimulation<-ListDataSimulation.MI2
		DataDesign<-DataDesign.MI2

	}
	
	
	if(n_data=="5000"){
		# ----- > loading simulated datasets
		load(paste0(Dir,"ListDataSimulation5000.MI2.RData"))
		load(paste0(Dir,"DataDesign5000.MI2.RData"))

		ListDataSimulation<-ListDataSimulation.MI2
		DataDesign<-DataDesign.MI2

	}
	
	
	if(n_data=="10000"){
		# ----- > loading simulated datasets
		load(paste0(Dir,"ListDataSimulation.MI2.RData"))
		load(paste0(Dir,"DataDesign.MI2.RData"))

		ListDataSimulation<-ListDataSimulation.MI2
		DataDesign<-DataDesign.MI2

	}
	
}


if(site=="Oeso"){
	
	# age
	load(file =paste0(Dir,"age.oeso.2000.Rdata" ))
	
	age.min<-round(quantile(age.oeso.2000,0.05)/5,0)*5
	cat(paste("\n", "\n", "minimum age","\n",sep=""))
	print(age.min)
	# ----- > loading formula for true model
	load(file =paste0(Dir,"modele.M0.oeso.Rdata" ))
	modele.true<-as.formula(modele.M0.oeso)
	cat(paste("\n", "\n", "true model","\n",sep=""))
	print(modele.true)
	
	form.true<-as.formula(~I(1/(1+temps)) + I(log(temps+1)) + temps + agex:(I(1/(1+temps)) + I(log(temps+1)) + temps) + agex + 
		I(agex^-0.5) + I(log(agex)))
		
	# ----- > loading true parameters
	load(file =paste0(Dir,"coef.M0.oeso.Rdata" ))
	coef.true<-coef.M0.oeso
	cat(paste("\n", "\n", "true coefs","\n",sep=""))
	print(coef.true)
	
	if(n_data=="2000"){
		# ----- > loading simulated datasets
		load(paste0(Dir,"ListDataSimulation.M0.RData"))
		load(paste0(Dir,"DataDesign.M0.RData"))

		ListDataSimulation<-ListDataSimulation.M0
		DataDesign<-DataDesign.M0

	}
	
	
	if(n_data=="5000"){
		# ----- > loading simulated datasets
		load(paste0(Dir,"ListDataSimulation5000.oeso.RData"))
		load(paste0(Dir,"DataDesign5000.oeso.RData"))

		ListDataSimulation<-ListDataSimulation.M0
		DataDesign<-DataDesign.M0

	}
	
	
	if(n_data=="10000"){
		# ----- > loading simulated datasets
		load(paste0(Dir,"ListDataSimulation.M0.RData"))
		load(paste0(Dir,"DataDesign.M0.RData"))

		ListDataSimulation<-ListDataSimulation.M0
		DataDesign<-DataDesign.M0

	}
	
}



#---------------------------------------------------------------------------------------
# FITTING MODELS AND CALCULATING INDICATORS
#---------------------------------------------------------------------------------------


cat(paste("\n", "\n", "Starting time, fit","\n",sep="")); print(date())
sink()

source(paste0(Dir,"00-Ajustement.R"))
sink(paste(Resultats,"/log.",site,".",n_data,".txt",sep=""),append=TRUE)
cat(paste("\n", "\n", "End time, fit","\n",sep="")); print(date())

cat("--------------------------------------------","\n")
cat("--------------------------------------------","\n")

cat(paste("\n", "\n", "Starting time, predictions","\n",sep="")); print(date())
sink()
source(paste0(Dir,"01-Calculs_Taux_Survie.R"))
sink(paste(Resultats,"/log.",site,".",n_data,".txt",sep=""),append=TRUE)
cat(paste("\n", "\n", "End time, predictions","\n",sep="")); print(date())

cat("--------------------------------------------","\n")
cat("--------------------------------------------","\n")

cat(paste("\n", "\n", "Starting time, RMISE","\n",sep="")); print(date())
sink()
source(paste0(Dir,"02-Calculs_Taux_MISE.R"))
sink(paste(Resultats,"/log.",site,".",n_data,".txt",sep=""),append=TRUE)
cat(paste("\n", "\n", "End time, RMISE","\n",sep="")); print(date())

cat("--------------------------------------------","\n")
cat("--------------------------------------------","\n")

heure.fin=Sys.time()
cat("\n", "\n")
print(heure.fin)

cat("\n", "\n", "Programme's duration","\n")
print(difftime(heure.fin, heure.debut,units = c("hours")))


sink()

# In the end you get dataframes containing every prediction needed to calculate the Bias, RMSE, RMISE and
# coverage probabilites.

# Concerning coverage probiblities you will need to eliminate the following files: 
# 821 for Oeso2000, 189 and 370 for Cervix2000 and 213 for Cervix5000. This needs to be done
# because, even if the model converged, the associated frequentist covariance matrices were not positive definite at convergence
		


#############################################################################


