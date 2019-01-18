\name{supGOGaU}
\alias{supGOGaU}

\title{Support of GOGaU Distribution}
\description{
This function provide useful support of GOGaU distribution for integration.
}

\usage{
supGOGaU(alpha, beta, a = 0, b = 1)
}
\arguments{
  \item{alpha}{The value of the first shape parameter. Must be positive and finite.}
  \item{beta}{The value of the second shape parameter. Must be positive and finite.}
  \item{a, b}{Lower and upper limits of the distribution. Must be finite.}
  }

\value{
  \item{lower}{The lower of GOGaU distribution useful support.}
  \item{upper}{The upper of GOGaU distribution useful support.}
 }
\references{
Bistoon Hosseini, Mahmoud Afshari, and Morad Alizadeh. "The Generalized Odd Gamma-G Family of Distributions:
Properties and Applications." Austrian Journal of Statistics 47.2 (2018): 69-89.
}
\author{Bistoon Hosseini, Mahmoud Afshari}
\examples{
supGOGaU(alpha = 0.5, beta = 1.5, a = 1, b = 4.5)
}

