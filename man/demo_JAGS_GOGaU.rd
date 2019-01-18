\name{demo_JAGS_GOGaU}
\alias{demo_JAGS_GOGaU}

\title{Demo of the JAGS GOGaU Module}
\description{
		This function provide a demo of the calculation of the Bayesian estimator using
		the GOGaU distribution module in the \pkg{JAGS} software. Here, 300 observations will generate from GOGaU(5, 0.28, 0, 10)
		and the Bayesian estimators are calculated by Gibbs sampling in \pkg{JAGS} software.
	}

\usage{
demo_JAGS_GOGaU()
}

\arguments{
This function has no parameter.
  }

\value{
	This function has no return value. The histogram of simulated data, The fitted GOGaU density of Bayesian estimators,
	and GOGaU(5, 0.28, 0, 10) density  is drawn.
 }
\references{
Bistoon Hosseini, Mahmoud Afshari, and Morad Alizadeh. "The Generalized Odd Gamma-G Family of Distributions:
Properties and Applications." Austrian Journal of Statistics 47.2 (2018): 69-89.
}
\author{Bistoon Hosseini, Mahmoud Afshari}

\examples{
demo_JAGS_GOGaU()
}

