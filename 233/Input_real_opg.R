##' ---
##' title: ""
##' date: ""
##' author: ""
##' output: pdf_document
##' ---


#+  message=FALSE
.libPaths("~/Rlibs/lib")
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(metrumrg)
library(parallel)
library(knitr)
opts_chunk$set(comment='.',echo=FALSE, message=FALSE)

mod <- mread("opg")
om <- omat(mod)
sg <- smat(mod)

mod %<>% mrgsolve:::collapse_omega() %>% mrgsolve:::collapse_sigma()

est <- c(TVCL=168,TVVC=2800,TVP1=443,TVP2=269,TVQ1=15.5,
         TVQ2=3.02,TVKA=0.0131,TVVMAX=13300,TVKM=6.74,TVFSC=0.0719)

rse <- c(3,2,16,14,16,13,4,13,11,0)
var <- (rse*est/100)^2

estt <- c(TVKSYN=0.864,TVKDEG=0.0204,TVIC50=5.38)
rsee <- c(8,6,21)
varr <- (rsee*estt/100)^2


theta <- c(est,estt)
Sigma <- diag(c(var,varr))
omega <- as.matrix(omat(mod))
sigma <- as.matrix(smat(mod))
dimnames(omega) <- list(NULL,NULL)
dimnames(sigma) <- list(NULL,NULL)

#+ echo=TRUE, comment='.'

##' # Fixed effect estimates:
#+ echo=FALSE
data_frame(Name = names(theta),
           Estimate = theta,
           variance = diag(Sigma)) %>% as.data.frame

##' # Between subject variability
om

##' # Residual error
sg

##' # Simulation of posterior
#+ echo=TRUE
simpost <- function(n) {
  post <- metrumrg::simpar(n,
                           theta,
                           covar=Sigma,
                           omega=omega,
                           sigma=sigma,
                           odf=100,sdf=1000)
  post <- post %>% as.data.frame
  nam <- names(post)
  nam <- sub(".", "", nam, fixed=TRUE)
  thetas <- which(grepl("TH",nam))
  nam[thetas] <- names(theta)
  nam <- sub("OM", "OMEGA",nam,fixed=TRUE)
  nam <- sub("SG", "SIGMA",nam,fixed=TRUE)
  names(post) <- nam
  post
}

#+ echo=TRUE
set.seed(111)
simpost(10)

#+ echo=TRUE
if(FALSE) {
  set.seed(22222)
  saveRDS(simpost(1000), file="opgpost.RDS")
}

##' # Session Info
sessionInfo()

