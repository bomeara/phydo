% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runanalysis.R
\name{render_phydo}
\alias{render_phydo}
\title{Create a file of results}
\usage{
render_phydo(taxon, format = "pdf", output_dir = getwd())
}
\arguments{
\item{taxon}{Clade of interest}

\item{format}{Format: pdf or html}

\item{output_dir}{Where to put the output; by default, current working directory}
}
\value{
Nothing, though a file is created in the current working directory
}
\description{
Create a file of results
}
\examples{
render_phydo("Tyto")
}
