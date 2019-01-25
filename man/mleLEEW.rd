\name{mleLEEW}
\alias{mleLEEW}

\title{The MLEs of LEEW}
\description{
This function calculates the MLEs of LEEW distribution.
}

\usage{
mleLEEW(x, model = "LEEW", fitplot = TRUE)
}
\arguments{
  \item{x}{The value of the first shape parameter. Must be finite.}
  \item{model}{"LEEW" (Log-Extended Exp-Weibull) or "LW" (Log-Weibull) model.}
  \item{fitplot}{Logical; if TRUE, histogram and fitted LEEW density is drawn.}
}

\value{
  \item{par}{The MLEs of parameters.}
  \item{loglike}{The log-likelihood value corresponding to MLEs.}
  \item{convergence}{An integer code. 0 indicates successful completion.
					To see possible error codes, see the \link{optim} function.
					}
  \item{KS}{Kolmogorov Smirnov test.}
  \item{AIC}{The value of Akaike Information Criterion.}
  \item{BIC}{The value of Bayesian Information Criterion.}
 }
\references{
Alizadeh, Morad, Mahmoud Afshari, Bistoon Hosseini, and Thiago G. Ramires.
"Extended exp-G family of distributions: Properties, applications and simulation."
Communications in Statistics-Simulation and Computation (2018): 1-16.
}
\author{Bistoon Hosseini, Mahmoud Afshari}
\examples{
set.seed(100)
x = runif(n = 200, 0, 10)
x = log(x)
mleLEEW(x, model = "LW", fitplot = TRUE)
mleLEEW(x, model = "LEEW", fitplot = TRUE)
}

