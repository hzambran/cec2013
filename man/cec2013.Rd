%% File cec2013.Rd
%% Part of the cec2013 R package, http://www.rforge.net/cec2013/ ; 
%%                                http://cran.r-project.org/web/packages/cec2013
%% Copyright 2013 Mauricio Zambrano-Bigiarini & Yasser Gonzalez Fernandez
%% Distributed under GPL 3 or later

\name{cec2013}
\alias{cec2013}
\title{Benchmark Functions for CEC 2013}
\description{
Common interface to the 28 CEC-2013 benchmark functions
}
\usage{
cec2013(i, x)
}

\arguments{
  \item{i}{numeric (integer) between 1 and 28, defining the number of the benchmark function to be evaluated on the \code{x} parameter set}
  \item{x}{Either a vector with 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 elements or a matrix with the same number of columns and one row for each vector to be evaluated}
}
\value{
numeric, with the value of the i-th CEC-2013 benchmark function evaluated in the vector \code{x} or for each row of \code{x} when the latter is a matrix
}
\references{
\cite{J. J. Liang, B-Y. Qu, P. N. Suganthan, Alfredo G. Hernandez-Diaz, \href{http://www.ntu.edu.sg/home/EPNSugan/index_files/CEC2013/Definitions\%20of\%20\%20CEC\%2013\%20benchmark\%20suite\%200117.pdf}{Problem Definitions and Evaluation Criteria for the CEC 2013 Special Session and Competition on Real-Parameter Optimization}, Technical Report 201212, Computational Intelligence Laboratory, Zhengzhou University, Zhengzhou China and  Technical Report, Nanyang Technological University, Singapore, January 2013}
}

\author{
Yasser Gonzalez Fernandez, \email{ygonzalezfernandez@gmail.com}, Mauricio Zambrano-Bigiarini, \email{mzb.devel@gmail.com}
}

\examples{
# CEC-2013 Benchmark function 1: Shifted Sphere, 10-dimensional space 
# x=[0,..0]^D:
func.num <- 1
D        <- 10
cec2013(func.num, rep(0, D))
}
