% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/molecule_serch.R
\name{molecule_serch}
\alias{molecule_serch}
\title{Screen molecules based on ingredients from Chinese herbal medicine}
\usage{
molecule_serch(molecule, ...)
}
\arguments{
\item{molecule}{The name of the molecule filter}

\item{...}{additional parameters}
}
\value{
data.frame
}
\description{
Screen molecules based on ingredients from Chinese herbal medicine
}
\examples{
molecule <- c("Hexyl acetate", "quercetin 3-o-rhamnopyranosyl")
molecule_filter<- molecule_serch(molecule)
head(molecule_filter)
}
