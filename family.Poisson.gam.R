#---------------------------------------------------------------------
# POISS.RS.SPLIT.R.GAM


# We adapt the Poisson family so it can handle expected mortality rates
POISS.RS.SPLIT.R.GAM <- function(link="log", MyData.dcatt)
{
	poissRS <- make.link("log")
	poissRS$linkfun <- function(mu, dcatt=MyData.dcatt)
	{
		if(any(mu<=dcatt))
		{
			print(summary(mu[mu<=dcatt]))
			print(summary(dcatt[mu<=dcatt]))
			cat(paste("\n",sum(mu<=dcatt)," / ",sum(mu<=dcatt)/length(mu),"\n",sep=""))
		}
		log(mu-dcatt)
	}
	poissRS$linkinv <- function(eta, dcatt=MyData.dcatt) exp(eta)+dcatt
	poissRS$mu.eta <- function(eta) exp(eta)

	# adding the second, third and fourth derivatives of the link function wrt the mean mu
	poissRS$d2link <- function(mu, dcatt=MyData.dcatt) -1/(mu-dcatt)^2
	poissRS$d3link <- function(mu, dcatt=MyData.dcatt) 2/(mu-dcatt)^3
	poissRS$d4link <- function(mu, dcatt=MyData.dcatt) -6/(mu-dcatt)^4
	
	# adding the first and second derivatives of the variance function
	poissRS$dvar <- function(mu) rep.int(1, length(mu))
	poissRS$d3var <- poissRS$d2var <- function(mu) rep.int(0, length(mu))

	# adding log saturated likelihood
	poissRS$ls <- function(y, w, n, scale)
	{
		res <- rep(0, 3)
		res[1] <- sum(dpois(y, y, log = TRUE) * w)
		res
	}

	poissRS$qf <- function(p, mu, wt, scale)
	{
		qpois(p, mu)
	}
	
	# adding qf and rd 
	poissRS$qf <- function(p, mu, wt, scale)
	{
		qpois(p, mu)
	}
	poissRS$rd <- function(mu, wt, scale)
	{
		rpois(length(mu), mu)
	}
 
	variance <- function(mu) mu
	validmu <- function(mu) all(mu > 0)
	dev.resids <- function(y, mu, wt) 2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
	aic <- function(y, n, mu, wt, dev) -2 * sum(dpois(y, mu, log = TRUE) * wt)

	initialize <- expression(
		{
			if (any(y < 0)) stop("negative values not allowed for the Poisson family")
			n <- rep.int(1, nobs)
			mustart <- y + 0.3
			
		})
	
	structure(list(family="poissonRS", link="poissonRS.link", linkfun=poissRS$linkfun, linkinv=poissRS$linkinv,
		d2link=poissRS$d2link, d3link=poissRS$d3link, d4link=poissRS$d4link, ls=poissRS$ls, 
		variance=variance, dvar=poissRS$dvar, d2var=poissRS$d2var, d3var=poissRS$d3var,
		dev.resids=dev.resids, aic=aic, qf=poissRS$qf, rd=poissRS$rd,
		mu.eta=poissRS$mu.eta, initialize=initialize, validmu=validmu,
		valideta=poissRS$valideta),class="family")
}



