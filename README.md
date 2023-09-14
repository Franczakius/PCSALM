# PCSALM
The R code to fit the parsimonious contaminated shifted asymmetric Laplace mixtures.

The files in this repository allow one to fit the parsimonious contaminated shifted asymmetric Laplace mixtures (PCSALM). The PCSALM is a family of mixture models that arise by imposing a combination of constraints on the scale matrix structure of a mixture of shifted asymmetric Laplace factor analyzers. 

The main function is PCSALM, which is available in PSCALM | AECM.R

The subfunctions called by PCSALM are available in PSCALM | Required Functions.R

Two illustrations of how to run the full family of models are given in PSCALM | Illustrations. R

There are three dependencies for this version of the code: library(Bessel), library(MixGHD), and library(timeR). If these packages are not installed, the user will receive an error before they get to the model fitting steps.

All files must be saved to the same folder on your machine, and the R working drive must be set to said folder. If this is not the case, an error will be generated indicating R can't find the required source files in PCSALM | Illustrations.R

The main function, PCSALM, has the following inputs:

x - a n by p matrix where the rows represent the observations and the columns represent the variables in the data set of interest

q - the number of latent factors 

G - the number of mixture components

init - the model parameters used to initialize the alternating-expectation conditional-maximization (AECM) algorithm (see PCSALM | Illustrations. R for examples)

max.it - the maximum number of iterations for the AECM algorithm used to fit the PCSALM

c.crit - the epsilon value for the convergence criteria used by the AECM algorithm

c.type - the type of convergence criteria, 1 is the difference between asymptotic estimate of the log-likelihood at iteration (k) and the log-likelihood value on iteration (k-1), 2 is the difference between asymptotic estimate of the log-likelihood 
at iteration (k) and the log-likelihood value on iteration (k), and 3 is the is the difference between asymptotic estimate of the log-likelihood at iteration (k) and iteration (k-1)

eta.cap - the maximum value allowed for eta, the degree of contamination

eta.C - a T or F value indicating whether or not eta should be constrained to be equal across groups

print.res - a T or F value indicating whether or not the results from the model fitting process should be printed
