
######################################################################################################
######################################################################################################
######################################################################################################

# This file contains code used in the Supplementary Material C.1 (only for sample size = 2,000)

	
	
# change this path (finish by /)
Dir <- "path/" 
setwd(Dir)



# warning : the bamlss models and the tensor_3 gss model have been commented out as they take a very long time to run. 
# Of course, you remain free of uncommenting them.


### packages
library(mgcv) # v1.8-26, we are using the gam function from mgcv. The bam function is generally faster but,
# as far as we know, it does not handle user-defined family as gam does so excess hazard models are not possible.
library(gss) # v2.1-9
library(rstpm2) # v1.4.5
library(survPen) # v1.0.1
library(R2BayesX) # v1.1-1
library(frailtypack) # v3.0.2.1
library(bamlss) # v1.0-1



# functions to split the data
source("fun.split.data.R")
source("Lexis.R")
source("family.Poisson.gam.R")

# expected rates (from the general French population)
rates <- dget(file="muaDF8933.dat"  )



# preparing data
data(datCancer) # from survPen package


#------------ split for mgcv

tempdata <- datCancer
tempdata$sexe <- 2 # cervix uteri cancer
tempdata$Ident <- 1:dim(tempdata)[1]


split2000 <- fun.split.data(tempdata, bands=c(seq(0,1,by=0.05),seq(1.1,5,by=0.1)), expected.rates=rates)
	
dim(split2000)	




# choose either 2 000, 20 000, or 100 000
num <- 2000



if (num==2000){
	don <- datCancer
	splitdon <- split2000
}


# age and year are centered
don$age <- trunc(don$age) - 70
don$yod <- trunc(don$yod) - 2000

splitdon$age <- trunc(splitdon$age) - 70
splitdon$yod <- trunc(splitdon$yod) - 2000
#----------------------


# for results

res.init <- rep(-1,8)

tab.summary <- data.frame(spline_time=res.init,tensor2=res.init,tensor3=res.init,tensor3_net=res.init)

rownames(tab.summary) <- c("frailtypack", "R2BayesX REML","R2BayesX MCMC", "bamlss", "gss", "rstpm2", "mgcv", "survPen")



sink(paste0(Dir,"compar_packages_",num,".lis"))
 
sessionInfo() 
 
 
cat("\n","dimensions of data, ", dim(don),"\n") 
 
cat("\n","dimensions of split data, ", dim(splitdon),"\n")  
 
#_____________________________________________________________________________


cat("\n","\n","Comparison of frailtypack, R2BayesX, bamlss, gss, rstpm2, mgcv and survPen", "\n")
cat("data = datCancer (from survPen package)", "\n")
cat("model = spline of time with df=10", "\n")
cat("type = overall survival", "\n")




# degrees of freedom
df1 <- 10


# bayesX REML
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","bayesX REML ","\n") 

time.bayesx <- system.time(fit.bayesx <- bayesx(dead~sx(fu, bs="bl",degree = 3, knots = df1), data=don,
family = "cox", method = "REML"))


print(summary(fit.bayesx))

cat("\n","bayex REML, execution time ", time.bayesx[3], "\n")

tab.summary[rownames(tab.summary)=="R2BayesX REML",c("spline_time")] <- time.bayesx[3]

rm(fit.bayesx)
gc()


# bayesX MCMC
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","bayesX MCMC ","\n") 

time.bayesx.MCMC <- system.time(fit.bayesx.MCMC <- bayesx(dead~sx(fu, bs="bl",degree = 3, knots = df1), data=don,
family = "cox", method = "MCMC"))


print(summary(fit.bayesx.MCMC))


cat("\n","bayex MCMC, execution time ", time.bayesx.MCMC[3], "\n")

tab.summary[rownames(tab.summary)=="R2BayesX MCMC",c("spline_time")] <- time.bayesx.MCMC[3]

rm(fit.bayesx.MCMC)
gc()

# bamlss
#cat("\n","\n","_____________________________________________________________________", "\n")
#cat("\n","bamlss ","\n") 

#f <- list(
#  Surv(fu, dead) ~ s(fu,bs="cr",k=df1)
  
#)

#set.seed(222)

#time.bamlss <- system.time(bam.fit <- bamlss(f, data = don, family = "cox",
#   subdivisions = 25, maxit = 1000,
#   n.iter = 6000, burnin = 3000, thin = 20, cores = 1,verbose=FALSE))


#print(summary(bam.fit))
	
	
#cat("\n","bamlss, execution time ", time.bamlss[3], "\n")

#tab.summary[rownames(tab.summary)=="bamlss",c("spline_time")] <- time.bamlss[3]

#rm(bam.fit)
#gc()
# frailtypack
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","frailtypack ","\n") 

time.frailtypack <- system.time(frailty.fit <- frailtyPenal(Surv(fu,dead)~1,data=don,n.knots=df1,kappa=10000,cross.validation=TRUE))



print(frailty.fit)

cat("\n","frailtypack, execution time ", time.frailtypack[3], "\n")

tab.summary[rownames(tab.summary)=="frailtypack",c("spline_time")] <- time.frailtypack[3]

rm(frailty.fit)
gc()
# gss
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","gss ","\n") 


time.gss <- system.time(gss.fit <- sshzd(Surv(fu,dead)~fu,data=don,nbasis=df1))


print(gss.fit)

cat("\n","gss, execution time ", time.gss[3], "\n")

tab.summary[rownames(tab.summary)=="gss",c("spline_time")] <- time.gss[3]

rm(gss.fit)
gc()
# survPen
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","survPen ","\n") 


time.survPen <- system.time(fit <- survPen(~smf(fu,df=df1),t1=fu,event=dead,data=don,method="LCV"))


print(summary(fit))

cat("\n","survPen, execution time ", time.survPen[3], "\n")

tab.summary[rownames(tab.summary)=="survPen",c("spline_time")] <- time.survPen[3]

rm(fit)
gc()
# rstpm2
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","rstpm2 ","\n") 


time.rstpm2 <- system.time(fit.rstpm2 <- pstpm2(Surv(fu,dead)~1, data=don, smooth.formula = ~s(log(fu),bs="cr",k=df1),
tvc = NULL,
bhazard = NULL,
sp=NULL,
criterion="GCV",
link.type="PH"))


print(summary(fit.rstpm2))

cat("\n","rstpm2, execution time ", time.rstpm2[3], "\n")

tab.summary[rownames(tab.summary)=="rstpm2",c("spline_time")] <- time.rstpm2[3]

rm(fit.rstpm2)
gc()
# mgcv
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","mgcv ","\n") 


time.mgcv <- system.time(mgcv.fit <- gam(dead~s(fu,bs="cr",k=df1),offset=log(tik),
		family=poisson,data=splitdon,scale=1,method="GCV.Cp" 
		))	  
			
		
print(summary(mgcv.fit))

cat("\n","mgcv, execution time ", time.mgcv[3], "\n")

tab.summary[rownames(tab.summary)=="mgcv",c("spline_time")] <- time.mgcv[3]


cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","the results so far ","\n") 
print(tab.summary)

rm(mgcv.fit)
gc()
#________________________________________________________________________________________________
#________________________________________________________________________________________________


cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")

#---- tensor

cat("\n","\n","Comparison of bamlss, gss, rstpm2, mgcv and survPen", "\n")
cat("data = datCancer (from survPen package)", "\n")
cat("model = tensor(time,age) df=c(5,5)", "\n")
cat("type = overall survival", "\n")




df2 <- c(5,5) 


gc()
# bamlss
#cat("\n","\n","_____________________________________________________________________", "\n")
#cat("\n","bamlss ","\n") 


#f <- list(
#  Surv(fu, dead) ~ te(fu,age,bs="cr",k=df2)
  
#)

#set.seed(222)

#time.bamlss2 <- system.time(bam.tensor <- bamlss(f, data = don, family = "cox",
#    subdivisions = 25, maxit = 1000,
#    n.iter = 6000, burnin = 3000, thin = 20, cores = 1,verbose=FALSE))

#print(summary(bam.tensor))

#cat("\n","bamlss, execution time ", time.bamlss2[3], "\n")


#tab.summary[rownames(tab.summary)=="bamlss",c("tensor2")] <- time.bamlss2[3]

#rm(bam.tensor)
#gc()
# gss
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","gss ","\n") 

time.gss2 <- system.time(gss.tensor <- sshzd(Surv(fu,dead)~fu*age,data=don,nbasis=prod(df2)))

print(gss.tensor)

cat("\n","gss, execution time ", time.gss2[3], "\n")

tab.summary[rownames(tab.summary)=="gss",c("tensor2")] <- time.gss2[3]

rm(gss.tensor)
gc()
# survPen
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","survPen ","\n") 

time.survPen2 <- system.time(fit.tensor <- survPen(~tensor(fu,age,df=df2),t1=fu,event=dead,data=don,method="LCV"))

print(summary(fit.tensor))

cat("\n","survPen, execution time ", time.survPen2[3], "\n")

tab.summary[rownames(tab.summary)=="survPen",c("tensor2")] <- time.survPen2[3]

rm(fit.tensor)
gc()
# rstpm2
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","rstpm2 ","\n") 

time.rstpm22 <- system.time(fit.rstpm2.tensor <- pstpm2(Surv(fu,dead)~1, data=don, smooth.formula = ~te(log(fu),age,bs="cr",k=df2),
tvc = NULL,
bhazard = NULL,
sp=NULL,
criterion="GCV",
link.type="PH"))

print(summary(fit.rstpm2.tensor))

cat("\n","rstpm2, execution time ", time.rstpm22[3], "\n")

tab.summary[rownames(tab.summary)=="rstpm2",c("tensor2")] <- time.rstpm22[3]

rm(fit.rstpm2.tensor)
gc()
# mgcv
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","mgcv ","\n") 


time.mgcv2 <- system.time(mgcv.tensor <- gam(dead~te(fu,age,bs="cr",k=df2),offset=log(tik),
		family=poisson,data=splitdon,scale=1,method="GCV.Cp" 
		))	  
			
		
print(summary(mgcv.tensor))

cat("\n","mgcv, execution time ", time.mgcv2[3], "\n")

tab.summary[rownames(tab.summary)=="mgcv",c("tensor2")] <- time.mgcv2[3]
rm(mgcv.tensor)
gc()
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","the results so far ","\n") 
print(tab.summary)


#________________________________________________________________________________________________
#________________________________________________________________________________________________


cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")


#---- tensor 3

cat("\n","\n","Comparison of bamlss, gss, rstpm2 and survPen", "\n")
cat("data = datCancer (from survPen package)", "\n")
cat("model = tensor(time,age,year) df=c(5,5,5)", "\n")
cat("type = overall survival", "\n")



df3 <- c(5,5,5)



# bamlss
#cat("\n","\n","_____________________________________________________________________", "\n")
#cat("\n","bamlss ","\n") 

#f <- list(
#  Surv(fu, dead) ~ te(fu,age,yod,bs="cr",k=df3)
  
#)

#set.seed(222)

#time.bamlss3 <- system.time(bam.tensor3 <- bamlss(f, data = don, family = "cox",
#    subdivisions = 25, maxit = 1000,
#    n.iter = 6000, burnin = 3000, thin = 20, cores = 1,verbose=FALSE))


#print(summary(bam.tensor3))
	
	
#cat("\n","bamlss, execution time ", time.bamlss3[3], "\n")


#tab.summary[rownames(tab.summary)=="bamlss",c("tensor3")] <- time.bamlss3[3]

#rm(bam.tensor3)
#gc()
# rstpm2
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","rstpm2 ","\n") 

time.rstpm23 <- system.time(fit.rstpm2.tensor3 <- pstpm2(Surv(fu,dead)~1, data=don, smooth.formula = ~te(log(fu),age,yod,bs="cr",k=df3),
tvc = NULL,
bhazard = NULL,
sp=NULL,
criterion="GCV",
link.type="PH"))

print(summary(fit.rstpm2.tensor3))

cat("\n","rstpm2, execution time ", time.rstpm23[3], "\n")

tab.summary[rownames(tab.summary)=="rstpm2",c("tensor3")] <- time.rstpm23[3]
rm(fit.rstpm2.tensor3)
gc()
# gss
#cat("\n","\n","_____________________________________________________________________", "\n")
#cat("\n","gss ","\n") 

#time.gss3 <- system.time(gss.tensor3 <- sshzd(Surv(fu,dead)~fu*age*yod,data=don,nbasis=prod(df3)))

#print(gss.tensor3)

#cat("\n","gss, execution time ", time.gss3[3], "\n")

#tab.summary[rownames(tab.summary)=="gss",c("tensor3")] <- time.gss3[3]
#rm(gss.tensor3)
#gc()
# survPen
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","survPen ","\n") 

time.survPen3 <- system.time(fit.tensor3 <- survPen(~tensor(fu,age,yod,df=df3),t1=fu,event=dead,data=don,method="LCV"))
time.survPen3


print(summary(fit.tensor3))

cat("\n","survPen, execution time ", time.survPen3[3], "\n")


tab.summary[rownames(tab.summary)=="survPen",c("tensor3")] <- time.survPen3[3]
rm(fit.tensor3)
gc()
# mgcv
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","mgcv ","\n") 


time.mgcv3 <- system.time(mgcv.tensor3 <- gam(dead~te(fu,age,yod,bs="cr",k=df3),offset=log(tik),
		family=poisson,data=splitdon,scale=1,method="GCV.Cp" 
		))	  
			
		
print(summary(mgcv.tensor3))

cat("\n","mgcv, execution time ", time.mgcv3[3], "\n")

tab.summary[rownames(tab.summary)=="mgcv",c("tensor3")] <- time.mgcv3[3]

cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","the results so far ","\n") 
print(tab.summary)
rm(mgcv.tensor3 )
gc()
#________________________________________________________________________________________________
#________________________________________________________________________________________________


cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")

# tensor 3 excess hazard

cat("\n","\n","Comparison of rstpm2, mgcv and survPen", "\n")
cat("data = datCancer (from survPen package)", "\n")
cat("model = tensor(time,age,year) df=c(5,5,5)", "\n")
cat("type = net survival", "\n")



df3 <- c(5,5,5)



# rstpm2
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","rstpm2 ","\n") 

time.rstpm24 <- system.time(fit.rstpm2.tensor3.net <- pstpm2(Surv(fu,dead)~1, data=don, smooth.formula = ~te(log(fu),age,yod,bs="cr",k=df3),
tvc = NULL,
bhazard = don$rate,
sp=NULL,
criterion="GCV",
link.type="PH"))

print(summary(fit.rstpm2.tensor3.net))

cat("\n","rstpm2, execution time ", time.rstpm24[3], "\n")

tab.summary[rownames(tab.summary)=="rstpm2",c("tensor3_net")] <- time.rstpm24[3]
rm(fit.rstpm2.tensor3.net)
gc()
# survPen

time.survPen4 <- system.time(fit.tensor3.net <- survPen(~tensor(fu,age,yod,df=df3),t1=fu,event=dead,data=don,method="LCV",expected=rate))
time.survPen4

print(summary(fit.tensor3.net))

cat("\n","survPen, execution time ", time.survPen4[3], "\n")

tab.summary[rownames(tab.summary)=="survPen",c("tensor3_net")] <- time.survPen4[3]
rm(fit.tensor3.net)
gc()
# mgcv
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","mgcv ","\n") 


time.mgcv4 <- system.time(mgcv.tensor3.net <- gam(dead~te(fu,age,yod,bs="cr",k=df3),offset=log(tik),
		family=POISS.RS.SPLIT.R.GAM(MyData.dcatt=splitdon$dcatt),data=splitdon,scale=1,method="GCV.Cp" 
		))	  

		
print(summary(mgcv.tensor3.net))

cat("\n","mgcv, execution time ", time.mgcv4[3], "\n")

tab.summary[rownames(tab.summary)=="mgcv",c("tensor3_net")] <- time.mgcv4[3]


cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")
cat("\n","\n","_____________________________________________________________________", "\n")

cat("\n","\n","Final results", "\n")

print(tab.summary)
rm(mgcv.tensor3.net)
gc()

sink()











