\name{LEEW}
\alias{dLEEW}
\alias{pLEEW}
\title{Log-Extended Exp-Weibull Distribution}
\description{
These functions provide density and distribution function of the Log-Extended Exp-Weibull 
distribution due to Alizadeh et al. (2018) specified by the pdf
\deqn{
f(y;\alpha,\beta,\mu,\sigma)=
\frac{\exp\left(\frac{y-\mu}{\sigma}-{\rm e}^{\frac{y-\mu}{\sigma}} \right)\left[1-\exp\left(-{\rm e}^{\frac{y-\mu}{\sigma}}\right) \right]^{\alpha-1}
\left\{\alpha +(\beta-\alpha)\left[1-\exp\left(-{\rm e}^{\frac{y-\mu}{\sigma}}\right)\right]^\beta \right\}}
{\sigma\left\{\left[1-\exp\left(-{\rm e}^{\frac{y-\mu}{\sigma}}\right)\right]^\alpha+1- \left[1-\exp\left(-{\rm e}^{\frac{y-\mu}{\sigma}}\right)\right]^\beta \right\}^2},
}
where \eqn{\alpha, \beta, \sigma > 0}.}
\usage{
dLEEW(x, alpha, beta, mu = 0, sigma = 1, log = FALSE)
pLEEW(x, alpha, beta, mu = 0, sigma = 1, log = FALSE)
}
\arguments{
  \item{x}{Scaler or vector of values at which the pdf or cdf needs to be computed}
  \item{alpha}{The value of the first shape parameter. Must be positive and finite.}
  \item{beta}{The value of the second shape parameter. Must be positive and finite.}
  \item{mu}{Value of mean. Must be finite.}
  \item{sigma}{Value of standard deviations. Must be positive and finite.}
  \item{log}{Logical; if TRUE, probabilities p are given as log(p).}
  }

\value{
An object of the same length as \code{x}, giving the pdf or cdf values computed at \code{x}.
 }
\references{
Alizadeh, Morad, Mahmoud Afshari, Bistoon Hosseini, and Thiago G. Ramires.
"Extended exp-G family of distributions: Properties, applications and simulation."
Communications in Statistics-Simulation and Computation (2018): 1-16.
}
\author{Bistoon Hosseini, Mahmoud Afshari}
\examples{
x = sort(rnorm(10, mean = -2 , sd = 1))
dLEEW(x, alpha = 0.5, beta = 1.2, mu = -2, sigma = 0.2, log = FALSE)
dLEEW(x, alpha = 0.5, beta = 1.2, mu = -2, sigma = 0.2, log = TRUE)
pLEEW(x, alpha = 0.5, beta = 1.2, mu = -2, sigma = 0.2)

