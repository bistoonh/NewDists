\name{supEEN}
\alias{supEEN}

\title{Support of EEN Distribution}
\description{
This function provide useful support of EEN distribution for integration.
}

\usage{
supEEN(alpha, beta, mu = 0, sigma = 1)
}
\arguments{
  \item{alpha}{The value of the first shape parameter. Must be positive and finite.}
  \item{beta}{The value of the second shape parameter. Must be positive and finite.}
  \item{mu}{Value of mean. Must be finite.}
  \item{sigma}{Value of standard deviations. Must be positive and finite.}
  }

\value{
  \item{lower}{The lower of EEN distribution useful support.}
  \item{upper}{The upper of EEN distribution useful support.}
 }
\references{
Alizadeh, Morad, Mahmoud Afshari, Bistoon Hosseini, and Thiago G. Ramires.
"Extended exp-G family of distributions: Properties, applications and simulation."
Communications in Statistics-Simulation and Computation (2018): 1-16.
}
\author{Bistoon Hosseini, Mahmoud Afshari}
\examples{
supEEN(alpha = 0.5, beta = 1.5, mu = 2, sigma = .5)
}

