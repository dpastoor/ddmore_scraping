##' ---
##' author: ""
##' date: ""
##' title: ""
##' output: 
##'   pdf_document:
##'     includes:  
##'       in_header: header.tex
##' ---


##' 
##' \head
##' 


#+ message=FALSE
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(knitr)
library(parallel)
opts_chunk$set(comment='.', message=FALSE)

##' Compile and load the OPG/NTX model
#+
mod <- mread("opg", end=528,delta=6)

##' A dosing event object for a single 3 mg/kg dose (assuming 70 kg individual)
#+
e <- ev(amt=70*3,cmt=1)

##' Execute the simulation
#+
sim <- 
  mod %>% Req(PKDV,NTX) %>% zero_re(sigma) %>%
  mrgsim(events=e,nid=50,atol=1E-30,obsonly=TRUE,add=0.01) %>%
  filter(time > 0)

#+
sim


#+ echo=FALSE
write.csv(file="Simulated_opg.txt", 
          as.data.frame(sim),
          quote=FALSE,row.names=FALSE)

##' 
##' \newpage
##'
##' ## Fc-OPG versus time
##' 
#+
ggplot(sim, aes(time,PKDV,group=ID)) + geom_line() + 
  scale_y_continuous(trans='log',breaks=10^seq(-4,4))

##' 
##' \newpage
##' 
##' ## uNTX versus time
##' 
#+
ggplot(sim, aes(time,NTX,group=ID)) + geom_line() 

##' 
##' \newpage
##'
##' # Simulation to evaluate 3 mg/kg dose 
##' 
##' - Simulate the median percent change from baseline at 2 weeks
##' after single 3 mg/kg dose
##' 

mod <- mread("opg")
mod %<>% mrgsolve:::collapse_omega() %>% mrgsolve:::collapse_sigma()

##' ### Read in our simulated posterior
#+
set.seed(770090)
post <- readRDS("opgpost.RDS") %>% sample_n(500)
omegas <- as_bmat(post, "OMEGA")
sigmas <- as_bmat(post,"SIGMA")

##' ### 3 mg/kg SC dose in N=200 70 kg patient
#+
sc3 <- expand.ev(amt=210,ID=1:200,IV=0)

#+
sim <- function(i,data,des) {
  mod %>%
    Req(PKDV,PDDV) %>%
    param(slice(post,i)) %>%
    omat(omegas[[i]]) %>% smat(sigmas[[i]]) %>%
    carry_out(dose,IV) %>%
    mrgsim(data=data,tgrid=des,obsonly=TRUE) %>%
    mutate(irep=i)
}

##' ### Simulate
##'   - If we're on a `unix-like` system, parallelize
##'   - If `windows`, just use regular `lapply`
#+
if(.Platform$OS.type=="windows") mclapply <- lapply

#+
mcRNG()
set.seed(2331991)
out <- mclapply(1:300, sim, data=sc3,
                des=tgrid(end=-1,add=c(0,336))) %>% bind_rows

##' ### Summarise the simulations:  
##'   - Group by `ID` and `irep`
##'   - Get the baseline `uNTX` observation 
##'   - Calculate percent change from baseline
##'   - Filter down to the week 2 observation
##'   - Summarize (median)
##' 
#+
sum <-
  out %>%
  group_by(ID,irep) %>%
  mutate(BASE = dplyr::first(PDDV),dDV = 100*(PDDV-BASE)/BASE) %>%
  ungroup %>%
  filter(time==336) %>%
  group_by(irep) %>%
  summarise(med = median(dDV))

##' ### The bottom line
##'   - The summarise the distribution of median percent change from baseline at week 2
#+
quantile(sum$med, c(0.025,0.5,0.975))

##'
##' \newpage
##'
#+
sessionInfo()



