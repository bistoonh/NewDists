dGOGaU = function(x,alpha, beta, a = 0, b = 1, log = FALSE)
	{
	par = c(alpha, beta, a, b)
	G = punif(x,par[3],par[4])
	g = dunif(x,par[3],par[4])
	Gb = G^par[2]
	d = par[2]*g*G^(par[1]*par[2]-1)*exp(-Gb/(1-Gb))/(gamma(par[1])*(1-Gb)^(par[1]+1))
	d[!is.finite(d)] = NA
	if(log == TRUE) d = log(d)
	return(d)
	}
	
dEEN = function(x, alpha, beta, mu = 0, sigma = 1, log = FALSE)
	{
	par = c(alpha, beta, mu, sigma)
	G = pnorm(x, par[3], par[4])
	g = dnorm(x, par[3], par[4])
	Ga = G^par[1]
	Gb = G^par[2]
	d = g * G^(par[1]-1) * (par[1] + (par[2]-par[1])*Gb) / (Ga + 1 - Gb)^2
	d[!is.finite(d)] = NA
	if(log == TRUE) d = log(d)
	return(d)
	}
	
dLEEW = function(x, alpha, beta, mu = 0, sigma = 1, log = FALSE)
	{
	z = (x - mu)/sigma
	u = 1 - exp(-exp(z))
	d = exp(z - exp(z)) * u^(alpha-1) * (alpha + (beta - alpha)*u^beta) / (sigma * (u^alpha + 1 - u^beta)^2)
	d[!is.finite(d)] = NA
	if(log == TRUE) d = log(d)
	return(d)
	}
	
pLEEW = function(x, alpha, beta, mu = 0, sigma = 1, log = FALSE)
	{
	z = (x - mu)/sigma
	u = 1 - exp(-exp(z))
	d = u^alpha / (u^alpha + 1 - u^beta)
	d[!is.finite(d)] = NA
	if(log == TRUE) d = log(d)
	return(d)
	}

pGOGaU = function(x, alpha, beta, a, b)
	{
	par = c(alpha, beta, a, b)
	Gb = punif(x, par[3], par[4])^par[2]
	p = pgamma(Gb/(1-Gb), par[1], 1)
	p[!is.finite(p)] = NA
	return(p)
	}

pEEN = function(x, alpha, beta, mu = 0, sigma = 1)
	{
	par = c(alpha, beta, mu, sigma)
	G = pnorm(x, par[3], par[4])
	g = dnorm(x, par[3], par[4])
	Ga = G^par[1]
	Gb = G^par[2]
	p = Ga / (Ga + 1 - Gb)
	p[!is.finite(p)] = NA
	return(p)
	}

qGOGaU = function(q, alpha, beta, a = 0, b = 1)
		{
		if(q<0 | 1<q) 
			{
			#print("q must be between [0,1].")
			return(NaN)
			} else{
				qs = qgamma(q, alpha, 1)
				qs = (qs/(1+qs))^(1/beta)
				qs = qunif(qs, a, b)
				return(qs)
				}
		}
qGOGaU = Vectorize(qGOGaU, "q")

qEEN = function(q, alpha, beta, mu = 0, sigma = 1)
	{
	if(q<0 | 1<q) 
		{
		#print("q must be between [0,1].")
		return(NaN)
		} else{
			f = function(t)
				{
				par = c(alpha, beta, mu, sigma)
				(q - pEEN(t, alpha, beta, mu, sigma))^2
				}
			qh = optim(mu, f, method = "L-BFGS-B", lower = mu-100*sigma, upper = mu+100*sigma)$par
			return(qh)
			}
	}
qEEN = Vectorize(qEEN, "q")

rGOGaU = function(n, alpha, beta, a = 0, b = 1)
	{
	par = c(alpha, beta, a, b)
	GI=rgamma(n,par[1],1)
	q = qunif((GI/(1+GI))^(1/par[2]),par[3],par[4])
	return(q)
	}
	
rEEN = function(n, alpha, beta, mu = 0, sigma = 1)
	{
	cfa = function(x, par) pGOGaU(x, par[1], par[2], par[3], par[4])
	par = c(alpha, beta, mu, sigma)
	q=runif(n, 0, 1)
	f <- function(P, fixed)
		{
		par <- fixed$par
		q <- fixed$q
		criterion <- q - cfa(P, par)
		criterion
		}
    P <- numeric(length(q))
    for(i in 1:length(q))
		{
		fixed <- list(par=par, q=q[i])
		root.p <- uniroot(f, lower = mu - 50*sigma, upper = mu + 50*sigma, fixed=fixed)
		P[i] <-root.p$root
		}
    P
	}

viewGOGaG = function()
	{
	fhGOGaG = function(par, xmax = 10, base = "Uniform")
		{
		if(base == "Uniform")
			{
			x <- seq(par[3], par[4], le = 10^3)
			G = punif(x,par[3],par[4])
			g = dunif(x,par[3],par[4])
			}
		if(base == "Weibull")
			{
			x <- seq(0, xmax, le = 10^3)
			G = pweibull(x,par[3],par[4])
			g = dweibull(x,par[3],par[4])
			}
		Gb = G^par[2]
		dd = par[2]*g*G^(par[1]*par[2]-1)*exp(-Gb/(1-Gb))/(gamma(par[1])*(1-Gb)^(par[1]+1))
		dd[!is.finite(dd)] = NA
		cd =  pgamma(Gb/(1-Gb),par[1],1)
		cd[!is.finite(cd)] = NA
		hd <- dd/(1 - cd)
		hd[!is.finite(hd)] = NA
		
		par(mfrow = c(1,2))
			plot(x, dd, type="l", lwd = 3, col = '#5B9BD5', ylab = "Density")
			if(base == "Uniform") plot(x, hd, type="l", lwd = 3, col = '#ED7D31', ylab = "Hazard", xlim = c(par[3], xmax))
			if(base == "Weibull") plot(x, hd, type="l", lwd = 3, col = '#ED7D31', ylab = "Hazard", xlim = c(0, xmax))
		par(mfrow = c(1,1))
		}
	manipulate(fhGOGaG(c(par1, par2, par3, par4), xmax, base),
				base = picker("Uniform", "Weibull", initial="Uniform"),
				xmax = slider(0, 30, initial = 1, label = "Max of x", step=.01),
				par1 = slider(.001, 5, initial = .5, label = "alpha", step=.01),
				par2 = slider(.001, 5, initial = .5, label = "beta", step=.01),
				par3 = slider(0, 10, initial = 0, label = "a", step=0.01),
				par4 = slider(0, 10, initial = 1, label = "b", step=0.01)
			)

	}

viewEEG = function()
	{
	fhEEG = function(par, xmax = 10, base = "Normal")
		{
		if(base == "Normal") 
			{
			x <- seq(par[3]-5*par[4], par[3]+5*par[4], le = 500)
			G = pnorm(x, par[3] ,par[4])
			g = dnorm(x,par[3],par[4])
			#G = pnorm(mpfr(x,200),par[3],par[4]); g = dnorm(x,par[3],par[4])
			}
		if(base == "Weibull")
			{
			x <- seq(0, xmax, le = 500)
			G = pweibull(x,par[3],par[4])
			g = dweibull(x,par[3],par[4])
			}
			
		pdf = function(par)
			{
			Ga = G^par[1]
			Gb = G^par[2]
			d = g * G^(par[1]-1) * (par[1] + (par[2]-par[1])*Gb) / (Ga + 1 - Gb)^2
			d[!is.finite(d)] = NA
			d
			}

		cdf = function(par)
			{
			Ga = G^par[1]
			Gb = G^par[2]
			d = Ga / (Ga + 1 - Gb)
			d[!is.finite(d)] = NA
			d
			}

		hrf <- function(par)
			{
			a <- pdf(par)
			b <- 1 - cdf(par)
			d <- a/b
			d[!is.finite(d)] = NA
			d
			}
		
		dd = pdf(par)
		hd = hrf(par)
	
		par(mfrow = c(1,2))
			plot(x, dd, type="l", lwd = 3, col = '#5B9BD5', ylab = "Density")
			plot(x, hd, type="l", lwd = 3, col = '#ED7D31', ylab = "Hazard", xlim = c(min(x), xmax))
		par(mfrow = c(1,1))
		}
	manipulate(fhEEG(c(par1, par2, par3, par4), xmax, base),
				base = picker("Normal", "Weibull", initial="Normal"),
				xmax = slider(0, 30, initial = 1, label = "Max of x", step=.01),
				par1 = slider(0, 10, initial = .5, label = "alpha", step=.01),
				par2 = slider(0, 10, initial = 1.5, label = "beta", step=.01),
				par3 = slider(-10, 10, initial = 0, label = "mu", step=0.01),
				par4 = slider(0, 10, initial = 1, label = "sigma", step=0.01)
			)

	}

viewLEEW = function()
	{
	fhLEEW = function(par)
		{
		x = seq(par[3] - 8 * par[4], par[3] + 5 *par[4], le = 1000)
		pdf = function(par1)
			{
			alpha = par[1]; beta = par[2]; mu = par[3]; sigma = par[4]
			z = (x - mu)/sigma
			u = 1 - exp(-exp(z))
			d = exp(z - exp(z)) * u^(alpha-1) * (alpha + (beta - alpha)*u^beta) / (sigma * (u^alpha + 1 - u^beta)^2)
			d[!is.finite(d)] = NA
			return(d)
			}
			
		cdf = function(par)
			{
			alpha = par[1]; beta = par[2]; mu = par[3]; sigma = par[4]
			z = (x - mu)/sigma
			u = 1 - exp(-exp(z))
			d = u^alpha / (u^alpha + 1 - u^beta)
			d[!is.finite(d)] = NA
			return(d)
			}


		S <- function(par)
			{
			d <- 1 - cdf(par)
			d[!is.finite(d)] = NA
			d
			}
		
		dd = pdf(par)
		Sd = S(par)
	
		par(mfrow = c(1,2))
			plot(x, dd, type="l", lwd = 3, col = '#5B9BD5', ylab = "Density")
			plot(x, Sd, type="l", lwd = 3, col = '#ED7D31', ylab = "Survival")
		par(mfrow = c(1,1))
		}
	manipulate(fhLEEW(c(par1, par2, par3, par4)),
				par1 = slider(0, 30, initial = 3, label = "alpha", step=.01),
				par2 = slider(0, 30, initial = 25, label = "beta", step=.01),
				par3 = slider(-10, 10, initial = 0, label = "mu", step=0.01),
				par4 = slider(0, 10, initial = 1, label = "sigma", step=0.01)
			)

	}

supGOGaU = function(alpha, beta, a = 0, b = 1)
		{
		m = matrix(NA, 10, 2)
		for(i in 1:10)
			{
			m[i, ] = qGOGaU(c(10^(-i), 1-10^(-i)), alpha, beta, a, b)
			}
		return(list(lower = min(m[, 1], na.rm = T), upper = max(m[, 2], na.rm = T)))
		}

supEEN = function(alpha, beta, mu = 0, sigma = 1)
		{
		m = matrix(NA, 10, 2)
		for(i in 1:10)
			{
			m[i, ] = qEEN(c(10^(-i), 1-10^(-i)), alpha, beta, mu, sigma)
			}
		return(list(lower = min(m[, 1], na.rm = T), upper = max(m[, 2], na.rm = T)))
		}
	
mGOGaU = function(alpha, beta, a = 0, b = 1)
	{
	subs = supGOGaU(alpha, beta, a = a, b = b)
	f <- function(x) x * dGOGaU(x, alpha, beta, a = a, b = b)
	meanf = integrate(f, subs$lower, subs$upper)$value
	return(meanf)		
	}

mEEN = function(alpha, beta, mu = 0, sigma = 1)
	{
	subs = supEEN(alpha, beta, mu = mu, sigma = sigma)
	f <- function(x) x * dEEN(x, alpha, beta, mu = mu, sigma = sigma)
	meanf = integrate(f, subs$lower, subs$upper)$value
	return(meanf)		
	}
	
vGOGaU = function(alpha, beta, a = 0, b = 1)
	{
	subs = supGOGaU(alpha, beta, a = a, b = b)
	f <- function(x) x^2 * dGOGaU(x, alpha, beta, a = a, b = b)
	varf = integrate(f, subs$lower, subs$upper)$value - mGOGaU(alpha, beta, a = a, b = b)^2
	return(varf)		
	}
	
vEEN = function(alpha, beta, mu = 0, sigma = 1)
	{
	subs = supEEN(alpha, beta, mu = mu, sigma = sigma)
	f <- function(x) x^2 * dEEN(x, alpha, beta, mu = mu, sigma = sigma)
	varf = integrate(f, subs$lower-.1*sigma, subs$upper+.1*sigma)$value - mEEN(alpha, beta, mu = mu, sigma = sigma)^2
	return(varf)		
	}
	
sGOGaU = function(alpha, beta, a = 0, b = 1)
	{
	subs = supGOGaU(alpha, beta, a = a, b = b)
	m1 = mGOGaU(alpha, beta, a = a, b = b)
	m2 = vGOGaU(alpha, beta, a = a, b = b)
	f <- function(x) (x-m1)^3 * dGOGaU(x, alpha, beta, a = a, b = b)
	sf = integrate(f, subs$lower, subs$upper)$value / m2^(3/2) 
	return(sf)
	}
	
sEEN = function(alpha, beta, mu = 0, sigma = 1)
	{
	subs = supEEN(alpha, beta, mu = mu, sigma = sigma)
	m1 = mEEN(alpha, beta, mu = mu, sigma = sigma)
	m2 = vEEN(alpha, beta, mu = mu, sigma = sigma)
	f <- function(x) (x-m1)^3 * dEEN(x, alpha, beta, mu = mu, sigma = sigma)
	sf = integrate(f, subs$lower-.15*sigma, subs$upper+.15*sigma)$value / m2^(3/2) 
	return(sf)
	}

kGOGaU = function(alpha, beta, a = 0, b = 1)
	{
	subs = supGOGaU(alpha, beta, a = a, b = b)
	m1 = mGOGaU(alpha, beta, a = a, b = b)
	m2 = vGOGaU(alpha, beta, a = a, b = b)
	f <- function(x) (x-m1)^4 * dGOGaU(x, alpha, beta, a = a, b = b)
	kf = integrate(f, subs$lower, subs$upper)$value / m2^2 
	return(kf)
	}
	
kEEN = function(alpha, beta,  mu = 0, sigma = 1)
	{
	subs = supEEN(alpha, beta, mu = mu, sigma = sigma)
	m1 = mEEN(alpha, beta, mu = mu, sigma = sigma)
	m2 = vEEN(alpha, beta, mu = mu, sigma = sigma)
	f <- function(x) (x-m1)^4 * dEEN(x, alpha, beta, mu = mu, sigma = sigma)
	kf = integrate(f, subs$lower-.2*sigma, subs$upper+.2*sigma)$value / m2^2 
	return(kf)
	}

entGOGaU = function(gamma, alpha, beta, a = 0, b = 1, explain = FALSE)
	{
	subs = supGOGaU(alpha, beta, a = a, b = b)
	if(gamma !=1)
		{
		f <- function(x) dGOGaU(x, alpha, beta, a = a, b = b)^gamma
		entf = log(integrate(f, subs$lower, subs$upper)$value)/(1-gamma)
		if(explain == TRUE) print(paste("The Rrnyi Entropy =", entf))
		} else{
				f <- function(x)
					{
					d = dGOGaU(x, alpha, beta, a = a, b = b)
					return(-log(d) * d)
					}
				entf = integrate(f, subs$lower, subs$upper)$value
				if(explain == TRUE) print(paste("The Shannon Entropy =", entf))
			  }
	return(entf)		
	}
entGOGaU = Vectorize(entGOGaU, "gamma")

entEEN = function(gamma, alpha, beta, mu = 0, sigma = 1, explain = FALSE)
	{
	subs = supEEN(alpha, beta, mu = mu, sigma = sigma)
	if(gamma !=1)
		{
		f <- function(x) dEEN(x, alpha, beta, mu = mu, sigma = sigma)^gamma
		entf = log(integrate(f, subs$lower-.1*sigma, subs$upper+.1*sigma)$value)/(1-gamma)
		if(explain == TRUE) print(paste("The Rrnyi Entropy =", entf))
		} else{
				f <- function(x)
					{
					d = dEEN(x, alpha, beta, mu = mu, sigma = sigma)
					return(-log(d) * d)
					}
				entf = integrate(f, subs$lower-.1*sigma, subs$upper+.1*sigma)$value
				if(explain == TRUE) print(paste("The Shannon Entropy =", entf))
			  }
	return(entf)		
	}
entEEN = Vectorize(entEEN, "gamma")

		
mleGOGaU = function(x, par0 = c(1, 1, min(x)-1, max(x)+1), a = NA, b = NA, fitplot = TRUE)
	{
	crGOGaU = function(dfa, cfa ,parh, x)
		{
		x_orderdenados = sort(x)
		v = cfa(x_orderdenados,parh)
		v[v == 1] = 1
		n = length(x)
		y = qnorm(v)
		u = pnorm((y - mean(y))/sqrt(var(y)))
		W_temp <- vector()
		A_temp <- vector()
		for (i in 1:n) {
			W_temp[i] = (u[i] - (2 * i - 1)/(2 * n))^2
			A_temp[i] = (2 * i - 1) * log(u[i]) + (2 * n + 1 - 
				2 * i) * log(1 - u[i])
		}
		A_2 = -n - mean(A_temp)
		W_2 = sum(W_temp) + 1/(12 * n)
		W_star = W_2 * (1 + 0.5/n)
		A_star = A_2 * (1 + 0.75/n + 2.25/n^2)
		p = length(parh)
		loglikefn = function(x, par)  -sum(log(dfa(x, par)))
		loglike = -1 * loglikefn(x, parh)
		AIC = -2 * loglike + 2 * p
		BIC = -2 * loglike + p * log(n)
		KS = ks.test(x = x, y = "cfa", par = as.vector(parh))
		return(list(W = W_star, A = A_star, KS = KS, AIC = AIC, BIC = BIC))
		}
	t = sort(x)
	if(is.na(a) & is.na(b))
		{
		lower = c(0.005, 0.005, -Inf, max(x)+.01)
		upper = c(Inf, Inf, min(x)-.01, Inf)
		loglike = function(pars) -sum(dGOGaU(x, pars[1], pars[2], pars[3], pars[4], log = TRUE))
		opt = optim(par0, loglike, method="L-BFGS-B", lower=lower, upper=upper)
		parh = opt$par
		if(fitplot == TRUE)
			{
			dd = dGOGaU(t, parh[1],  parh[2], parh[3], parh[4])
			hist(x, prob = T, ylim=c(0,1.1*max(dd, na.rm = T)))
			lines(t, dd, lwd = 2, col = '#5B9BD5')
			}
		dfa = function(x, parh) dGOGaU(x, parh[1], parh[2], parh[3], parh[4])
		cfa = function(x, parh) pGOGaU(x, parh[1], parh[2], parh[3], parh[4])
		crs = crGOGaU(dfa=dfa, cfa=cfa, parh=parh, x = x)
		} else if(!is.na(a) & is.na(b)){
					if(3 < length(par0))par0 = par0[-3]
				 	lower = c(0.01, 0.01, max(x)+.01)
					upper = c(Inf, Inf, Inf)
					loglike = function(pars) -sum(dGOGaU(x, pars[1], pars[2], a, pars[3], log = TRUE))
					opt = optim(par0, loglike, method="L-BFGS-B", lower=lower, upper=upper)
					parh = opt$par
					if(fitplot == TRUE)
						{
						dd = dGOGaU(t, parh[1],  parh[2], a, parh[3])
						hist(x, prob = T, ylim = c(0,1.1*max(dd, na.rm = T)))
						lines(t, dd, lwd = 2, col = '#5B9BD5')
						}
					dfa = function(x, parh) dGOGaU(x, parh[1],  parh[2], a, parh[3])
					cfa = function(x, parh) pGOGaU(x, parh[1],  parh[2], a, parh[3])
					crs = crGOGaU(dfa=dfa, cfa=cfa, parh=parh, x = x)
			} else if(is.na(a) & !is.na(b)){
					if(3 < length(par0))par0 = par0[-4]
					lower = c(0.01, 0.01, -Inf)
					upper = c(Inf, Inf, min(x)-.001)
					loglike = function(pars) -sum(dGOGaU(x, pars[1], pars[2], pars[3], b, log = TRUE))
					opt = optim(par0, loglike, method="L-BFGS-B", lower=lower, upper=upper)
					parh = opt$par
					if(fitplot == TRUE)
						{
						dd = dGOGaU(t, parh[1],  parh[2], parh[3], b)
						hist(x, prob = T, ylim = c(0,1.1*max(dd, na.rm = T)))
						lines(t, dd, lwd = 2, col = '#5B9BD5')
						}
					dfa = function(x, parh) dGOGaU(x, parh[1],  parh[2], parh[3], b)
					cfa = function(x, parh) pGOGaU(x, parh[1],  parh[2], parh[3], b)
					crs = crGOGaU(dfa=dfa, cfa=cfa, parh=parh, x = x)
				} else {
						if(2 < length(par0)) par0 = par0[-(3:4)]
						lower = c(0.01, 0.01)
						upper = c(Inf, Inf)
						loglike = function(pars) -sum(dGOGaU(x, pars[1], pars[2], a, b, log = TRUE))
						opt = optim(par0, loglike, method="L-BFGS-B", lower=lower, upper=upper)
						parh = opt$par
						if(fitplot == TRUE)
							{
							dd = dGOGaU(t, parh[1],  parh[2], a, b)
							hist(x, prob = T, ylim = c(0,1.1*max(dd, na.rm = T)))
							lines(t, dd, lwd = 2, col = '#5B9BD5')
							}
						dfa = function(x, parh) dGOGaU(x, parh[1],  parh[2], a, b)
						cfa = function(x, parh) pGOGaU(x, parh[1],  parh[2], a, b)
						crs = crGOGaU(dfa=dfa, cfa=cfa, parh=parh, x = x)
						}
	
		return(list(par = opt$par, loglike = -opt$value, convergence = opt$convergence, 
					W = crs$W, A = crs$A, KS = crs$KS, AIC = crs$AIC,  BIC = crs$BIC))
		}

mleEEN = function(x, par0 = c(1, 1, mean(x), sd(x)), fitplot = TRUE)
	{
	crEEN = function(dfa, cfa ,parh, x)
		{
		x_orderdenados = sort(x)
		v = cfa(x_orderdenados,parh)
		v[v == 1] = 1
		n = length(x)
		y = qnorm(v)
		u = pnorm((y - mean(y))/sqrt(var(y)))
		W_temp <- vector()
		A_temp <- vector()
		for (i in 1:n) {
			W_temp[i] = (u[i] - (2 * i - 1)/(2 * n))^2
			A_temp[i] = (2 * i - 1) * log(u[i]) + (2 * n + 1 - 
				2 * i) * log(1 - u[i])
		}
		A_2 = -n - mean(A_temp)
		W_2 = sum(W_temp) + 1/(12 * n)
		W_star = W_2 * (1 + 0.5/n)
		A_star = A_2 * (1 + 0.75/n + 2.25/n^2)
		p = length(parh)
		loglikefn = function(x, par)  -sum(log(dfa(x, par)))
		loglike = -1 * loglikefn(x, parh)
		AIC = -2 * loglike + 2 * p
		BIC = -2 * loglike + p * log(n)
		KS = ks.test(x = x, y = "cfa", par = as.vector(parh))
		return(list(W = W_star, A = A_star, KS = KS, AIC = AIC, BIC = BIC))
		}
	t = sort(x)
	loglike = function(pars) -sum(dEEN(x, pars[1], pars[2], pars[3], pars[4], log = TRUE))
	opt = optim(par0, loglike, method = 'Nelder-Mead', control = list(maxit = 5000))
	parh = opt$par
	if(fitplot == TRUE)
		{
		dd = dEEN(t, parh[1],  parh[2], parh[3], parh[4])
		hist(x, prob = T, ylim=c(0,1.1*max(dd, na.rm = T)))
		lines(t, dd, lwd = 2, col = '#5B9BD5')
		}
	dfa = function(x, parh) dEEN(x, parh[1], parh[2], parh[3], parh[4])
	cfa = function(x, parh) pEEN(x, parh[1], parh[2], parh[3], parh[4])
	crs = crEEN(dfa=dfa, cfa=cfa, parh=parh, x = x)
	return(list(par = opt$par, loglike = -opt$value, convergence = opt$convergence, 
					W = crs$W, A = crs$A, KS = crs$KS, AIC = crs$AIC,  BIC = crs$BIC))
		}

		
mleLEEW = function(x, model = "LEEW", fitplot = TRUE)
	{
	crLEEW = function(dfa, cfa ,parh, x)
		{
		x_orderdenados = sort(x)
		v = cfa(x_orderdenados,parh)
		v[v == 1] = 1
		n = length(x)
		p = length(parh)
		loglikefn = function(x, par)  -sum(log(dfa(x, par)))
		loglike = -1 * loglikefn(x, parh)
		AIC = -2 * loglike + 2 * p
		BIC = -2 * loglike + p * log(n)
		KS = ks.test(x = x, y = "cfa", par = as.vector(parh))
		return(list(KS = KS, AIC = AIC, BIC = BIC))
		}
	if(model == "LEEW") 
		{
		logdfa = function(y, par)
			{
			alpha = par[1]; beta = par[2]; mu = par[3]; sigma = par[4]
			z = (y - mu)/sigma
			u1 = 1 - exp(-exp(z))
			ret = (z - exp(z) + (alpha - 1) * log(u1) + log(alpha + (beta - alpha) * u1^beta) -
					log(sigma)- 2 * log(u1^alpha + 1 - u1^beta))
			ret[!is.finite(ret)] = NA
			return(ret)
			}
		findinit = function(r, y)
			{
			par0 = rep(NA, 4)
			s = 0
			i = 1
			while(s == 0 & i < r)
				{
				par0 = c(runif(2, 0, 100), runif(1, mean(y) - 1*sd(y), mean(y) - 1*sd(y)), runif(1, .5*sd(y), 2*sd(y)))
				if(!is.na(logdfa(y, par0))) s = 1 else i = i+1
				}
			if(i == r) return(rep(NA, 4)) else return(par0)
			}
		bestopt = function(r = 100, y)
			{
			pars = t(replicate(r, findinit(10^3, y)))
			f = function(par) -sum(logdfa(y, par))
			optimvals = function(par0) optim(par0, f, method = "Nelder-Mead")$value
			values = apply(pars, 1, optimvals)
			t = which.min(values)
			opt = optim(pars[t, ], f, method = "Nelder-Mead", control = list(maxit = 5000), hessian = T)
			sdpar = sqrt(diag(-solve(-opt$hessian)))
			return(list(par = opt$par, value = opt$value, convergence = opt$convergence, sdpar = sdpar))
			}
		opt = bestopt(r = 50, x)
		parh = opt$par
		if(fitplot == TRUE)
			{
			t = sort(x)
			dd = exp(logdfa(t, c(parh[1],  parh[2], parh[3], parh[4])))
			hist(x, prob = T, ylim=c(0,1.1*max(dd, na.rm = T)))
			lines(t, dd, lwd = 2, col = '#5B9BD5')
			}
		dfa = function(x, parh) dLEEW(x, parh[1], parh[2], parh[3], parh[4])
		cfa = function(x, parh) pLEEW(x, parh[1], parh[2], parh[3], parh[4])
		crs = crLEEW(dfa=dfa, cfa=cfa, parh=parh, x = x)	
		}
	if(model == "LW")	
		{
		logdfa = function(y, par)
			{
			alpha =1; beta = 1; mu = par[1]; sigma = par[2]
			z = (y - mu)/sigma
			u1 = 1 - exp(-exp(z))
			ret = (z - exp(z) + (alpha - 1) * log(u1) + log(alpha + (beta - alpha) * u1^beta) -
					log(sigma)- 2 * log(u1^alpha + 1 - u1^beta))
			ret[!is.finite(ret)] = NA
			return(ret)
			}
		findinit = function(r, y)
			{
			par0 = rep(NA, 2)
			s = 0
			i = 1
			while(s == 0 & i < r)
				{
				par0 = c(runif(1, mean(y) - 1*sd(y), mean(y) - 1*sd(y)), runif(1, .5*sd(y), 2*sd(y)))
				if(!is.na(logdfa(y, par0))) s = 1 else i = i+1
				}
			if(i == r) return(rep(NA, 2)) else return(par0)
			}
		bestopt = function(r = 100, y)
			{
			pars = t(replicate(r, findinit(10^3, y)))
			f = function(par) -sum(logdfa(y, par))
			optimvals = function(par0) optim(par0, f, method = "Nelder-Mead", control=list(maxit=5000))$value
			values = apply(pars, 1, optimvals)
			t = which.min(values)
			opt = optim(pars[t, ], f, method = "Nelder-Mead", control = list(maxit = 5000), hessian = T)
			sdpar = sqrt(diag(-solve(-opt$hessian)))
			return(list(par = opt$par, value = opt$value, convergence = opt$convergence, sdpar = sdpar))
			}
		opt = bestopt(r = 20, x)
		parh = opt$par
		if(fitplot == TRUE)
			{
			t = sort(x)
			dd = exp(logdfa(t, c(parh[1],  parh[2])))
			hist(x, prob = T, ylim=c(0,1.1*max(dd, na.rm = T)))
			lines(t, dd, lwd = 2, col = '#5B9BD5')
			}
		dfa = function(x, parh) dLEEW(x, 1, 1, parh[1], parh[2])
		cfa = function(x, parh) pLEEW(x, 1, 1, parh[1], parh[2])
		crs = crLEEW(dfa=dfa, cfa=cfa, parh=parh, x = x)	
		}

	return(list(par = opt$par, loglike = -opt$value, convergence = opt$convergence, 
					KS = crs$KS, AIC = crs$AIC,  BIC = crs$BIC))
	}
		
demo_JAGS_GOGaU = function()
	{
	load.module("GOGaU")

	# JAGS model
	model.string =  
		"
		model
			{
			alpha ~ dgamma(0.001, 0.001)
			beta ~ dgamma(0.001, 0.001)
			a ~ dunif(-40,minx)
			b ~ dunif(maxx,50)

			for (i in 1:N)
				{
				x[i] ~ dGOGaU(alpha, beta, a, b)
				}
			}
		"

	# simulation data
	coli = c('#5B9BD5', '#ED7D31', '#9B9B9B', '#70AD47', '#FFC000', '#A532A5', '#4076C4', '#255E91')
	N = 300
	set.seed(100)
	par = c(5, 0.28, 0, 10)
	x <- rGOGaU(N, par[1], par[2], par[3], par[4])

	dat <- list(x=x, N=N, minx = min(x), maxx = max(x))

	# inits
	inits1 <- list(alpha = par[1] , beta = par[2], a=min(x)-1, b = max(x)+1)
	inits <- list(inits1)

	# sample
	model.spec <- textConnection(model.string)
	j.model = NA
	j.model <- jags.model(model.spec, data = dat, inits=inits, n.chains=1, n.adapt=1000)
	j.samples <- coda.samples(j.model, c("alpha", "beta", "a", "b"), n.iter=6000, thin=1)


	# plot
	s = summary(j.samples)
	jagsparh = s$statistics[c("alpha", "beta", "a", "b"), 1]
	d1 = dGOGaU(sort(x), par[1], par[2], par[3], par[4])
	d2 = dGOGaU(sort(x), jagsparh[1], jagsparh[2], jagsparh[3], jagsparh[4])
	hist(x, prob = T, ylim = c(0, max(d1,d2)))
		lines(sort(x), d1, lwd = 2, col = '#5B9BD5')
		lines(sort(x), d2, lwd = 2, col = '#ED7D31')
		legend(0.5, max(d1,d2),
				c('Real density', 'Bayes estimator'),
				col = c('#5B9BD5', '#ED7D31'),
				lwd = c(2, 2)
				) 
	}
	
