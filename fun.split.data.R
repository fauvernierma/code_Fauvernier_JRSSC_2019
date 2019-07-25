#---------------------------------------------------------------------
# fun.split.data

# We specify a dataframe and a list of bands of follow-up times.
# From these bands, the function creates a split version of the original dataframe
fun.split.data <- function(dat, bands, expected.rates)
{

# preparing data with lexis function
splitdata <- lexis(entry = 0, exit = fu, fail = dead, breaks = bands,
                         include = list(Ident,age,yod,sexe),
                          data = dat)

# length of each time interval						  
splitdata$tik      <-  splitdata$Exit - splitdata$Entry

# 'point milieu' approximation of the cumulative hazard by taking the middle of each interval as
# the observed follow-up time
splitdata$fu <- (splitdata$Entry + splitdata$Exit) / 2

# redefining the name of the event indicator
splitdata$dead <- splitdata$Fail

# if death does occur, the final observed time must correspond to the true time of death and not
# the middle of the final interval
splitdata[splitdata$dead == 1, c("fu")] <- splitdata[splitdata$dead == 1, c("Exit")]
 
 
# merging with expected rates table using age and year of death
splitdata$age.current <- floor(floor(splitdata$age) + splitdata$fu)
splitdata$yod.current <- floor(floor(splitdata$yod) + splitdata$fu)


splitdata <- merge(splitdata,expected.rates,by.x=c("yod.current","age.current","sexe"),by.y=c("ANNEE","AGEX","SEXE"))

# expected numbers of deaths
splitdata$dcatt=splitdata$tik*splitdata$MUA  
 
 
return(splitdata)

}
