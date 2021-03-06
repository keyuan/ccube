% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccube.R
\name{MergeClusters}
\alias{MergeClusters}
\title{Merge clusters}
\usage{
MergeClusters(res = res, ssm = ssm, tol = 1e-08, maxiter = 100,
  epi = 0.001)
}
\arguments{
\item{res}{ccube results list}

\item{ssm}{ccube data frame}

\item{tol}{stopping condition in VBEM}

\item{maxiter}{maximum iteration in VBEM}

\item{epi}{sequencing error}
}
\value{
res ccube results list
}
\description{
Merge clusters
}
