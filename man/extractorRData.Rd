% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{extractorRData}
\alias{extractorRData}
\title{Extract an object from a .RData file}
\usage{
extractorRData(file, object)
}
\arguments{
\item{file}{A character string specifying the path to the RData file.}

\item{object}{A character string specifying the name of the object to extract.}
}
\value{
The object extracted from the RData file.
}
\description{
This function extracts an object from an RData file created by R's save() command.
}
