############################################################################################################
### A COLLECTION OF FUNCTIONS FOR USE IN SURIVIAL / TIME TO EVENT ANALYSIS, SIMULATION AND VISUALISATION ###
############################################################################################################

library(deSolve)
library(survival)
library(plyr)

############################################################################################################
### FUNCTION LIST ##########################################################################################
############################################################################################################ 
#   Models_summary()                              creates a table of modelling results and basic GOF plots
#   Hazard()                                      model parameters are supplied to the function and the hazard function is returned
#   CumHaz()                                      hazard function supplied by Hazard() is integrated to return the cumulative hazard
#   IntHaz()                                      a convenient wrapper for Hazard() and CumHaz() for calculating survival curve from a set of model parameters
#   CompeteEvents()                               carries out analysis of competing events leading to an ajustment of simulated events to give more realistic results in presence of competing events
#   calcEventVPC()                                carries out calculations for production of an event VPC
#   plotEventVPC()                                for plotting of event VPC
#   plot_KM()                                     for generation of Kaplan-Meier plots from data frame of event times and event types                                
#   add_KM_conf_bands()                           for adding confidence bands from multiple simulations to an existing KM plot
#   sim_events()                                  for simulation of an event for each patient in a study
#   repeated_sim_events()                         function to produce replicate simulations of events for each patient in a study
#   calc_KM_conf_bands()                          function to calculate confidence interval as central X% range of many simulated Kaplan-Meier curves
#   adjust_events_to_appointment_times()          function to adjust simulated event times to the next schecduled observation appointment for each patient
#   adjust_repeat_events_to_appointment_times()   function to adjust repeated simulated event times to the next scheduled observation appointment for each patient
############################################################################################################
############################################################################################################


##############################################################
### FUNCTION FOR COMPARING MODELS ############################
##############################################################
### LIST OF FUNCTION ARGUMENTS
#
# folders                    # vector of folders that contain the finished NONMEM runs of the models of interest.
# model_names                # vector of model names
# ref_model                  # reference model index of "folders" for calculating DOFV and AIC
# col_surv                   # column of .tab NONMEM file that contains values of the survival curve
# outputfolder               # output folder the .csv and .png files will be saved to
# create_plot                # if TRUE, a KM plot together with survival curves from all model outputs are plotted
#
##############################################################

Models_summary = function( folders, model_names, ref_model, col_surv, outputfolder, output_file_name, create_plot=FALSE, ... ){
  
  #########################################
  # SUMMARY FUNCTIONS
  #########################################
  get_DOFV_DAIC = function( BaseOFV, BaseThetas, CurrentOFV, CurrentThetas){
    
    DOFV = CurrentOFV - BaseOFV
    
    BaseAIC = BaseOFV + 2*BaseThetas
    CurrentAIC = CurrentOFV + 2*CurrentThetas  
    
    DAIC = CurrentAIC - BaseAIC
    
    return( c(DOFV, DAIC) )
  }
  
  # function to check if minimisation was successful based
  # on the res.csv file
  checkMinRes = function( res_dat ){
    
    #
    if (any(grepl("Parameter", colnames(res_dat)))==FALSE){
      return(0)
    }
    # rows of theta
    ths = grep("THETA(", res_dat[,"Parameter"], fixed=TRUE)
    if (length(ths)==0){
      return(0)
    }
    
    return(1)
  }
  
  # function for model summary
  model_summary = function( folder, model, base_summary=NULL ){
    
    model_summary = data.frame( RunRef = tail(unlist(strsplit(folder, "/", fixed=TRUE)),1) )
    
    model_summary$Model = model
    
    model_summary$Converged = 0
    if (any(grepl("res.csv", list.files(folder)))==FALSE){
      warning("No result file in ", folder, "\n")
      return(model_summary)
    }
    res_dat = read.csv(paste0(folder,"res.csv"))  
    
    if ( checkMinRes(res_dat) == 0 ){
      warning( "Minimisation not successfull ", folder, "\n")
      return(model_summary)
    }
    
    
    
    model_summary$Hit_bound = ifelse(length(grep(pattern = "zero gradients found", res_dat[,1], ignore.case = TRUE))>0, 1, 0)
    model_summary$Converged = ifelse(length(grep(pattern = "0MINIMIZATION SUCCESSFUL", res_dat[,1], ignore.case = TRUE))>0, 1, 0)
    
    OFV = as.numeric(as.character(res_dat[which(res_dat$Parameter=="OFV")[1],"Estimate"]))
    model_summary$OFV = OFV
    
    # all theta information
    ths = grep("THETA(", res_dat[,"Parameter"], fixed=TRUE)
    ths_len = length(ths)
    
    th_nms  = as.character((res_dat[ths,"Parameter"]))
    th_vls  = as.numeric(as.character((res_dat[ths,"Estimate"])))
    th_ses  = as.numeric(as.character((res_dat[ths,"Error"])))
    
    for (j in 1:ths_len){
      th_nms[j] = paste0("THETA",j)
    }
    
    if (length(base_summary)==0){
      BaseOFV  = OFV
      BaseThetas = length(ths)
    } else {
      BaseOFV  = base_summary$OFV
      BaseThetas = base_summary$n
    }
    
    DOFV_DAIC = get_DOFV_DAIC( BaseOFV, BaseThetas, OFV, length(ths) )
    
    model_summary$n    = ths_len
    model_summary$DOFV = DOFV_DAIC[1]
    model_summary$AIC = DOFV_DAIC[2]
    
    for ( j in 1:ths_len){
      model_summary[ , th_nms[j] ] = th_vls[j]
      
      model_summary[ , paste0(th_nms[j],"_SE") ] = th_ses[j]
      
      if (is.na(th_ses[j]) ){
        model_summary[ , paste0(th_nms[j],"_95%") ] = NA
      } else {
        temp = c(th_vls[j] - 1.96*th_ses[j], th_vls[j] + 1.96*th_ses[j])
        model_summary[ , paste0(th_nms[j],"_95%") ] = paste0("(", signif(temp[1],3), " ; ", signif(temp[2],3), ")")
      }
    }
    
    return(model_summary)
    
  }
  #########################################
  #########################################
  
  
  #########################################
  # COMMON ERROR TYPES
  #########################################
  if ( length(folders) != length(model_names) ) stop("Number of folders and model names do not match")
  
  if ( !is.numeric(ref_model) | ref_model<=0 ) stop("ref_model must be numeric and >0")
  
  if ( ref_model>length(folders) ) stop("ref_model should be less than or equal to the number of different models")
  #########################################
  #########################################
  
  # Rearranging such that reference model is always the first
  folders = replace(folders, c(1, ref_model), folders[c(ref_model, 1)])
  model_names = replace(model_names, c(1, ref_model), model_names[c(ref_model, 1)])
  
  
  #########################################
  # GET AIC AND TAB FILES
  #########################################
  # go through each folder get AIC and extract tab file
  for (i in 1:length(folders)){
    
    # gets the summary
    if (i==1){
      all_model_summary = model_summary(folders[i], model_names[i], NULL)
      if ( all_model_summary$Converged == 0 )  stop("Reference model not converged")
      model_names[i] = paste0(model_names[i], "; AIC: ", round(tail(all_model_summary$AIC,1),3) )
    } else {
      new_model_summary = model_summary(folders[i], model_names[i],all_model_summary[1,])
      all_model_summary = rbind.fill(all_model_summary,new_model_summary)
      model_names[i] = paste0(model_names[i], "; AIC: ", round(tail(all_model_summary$AIC,1),3 ) )
    }
    
    # checks if there is a tab file in foler
    if (any(grepl(".tab", list.files(folders[i])))==FALSE){
      warning("No tab file in ", folders[i], "\n")
      next
    }
    # get tabfile name
    tabfile = list.files(folders[i])[grepl(".tab", list.files(folders[i]))]
    
    # Read in tab file. If more than 1 .tab then read only the first one
    d = read.table(paste0(folders[i],tabfile[1]),header=TRUE,sep="",stringsAsFactors=FALSE,skip=1)
    d = d[ order(d$TIME),]
    
    # base model for KM plot is stored first
    if ( i==1 ){
      all_tab = list()
      all_tab[[i]] = d[ d$EVID==0, c("TIME", "DV")]
    }
    
    all_tab[[i+1]] = d[ d$EVID==0, c("TIME", col_surv)]
    
    
  }
  print(all_model_summary)
  try(write.csv(all_model_summary,paste0(outputfolder, output_file_name,".csv",sep=""),quote = FALSE,row.names = FALSE))
  #########################################
  #########################################
  
  
  
  #########################################
  # PLOTTING
  #########################################
  
  if (create_plot==TRUE){
    cols=c("blue","red","green","purple", "violet", "gold","brown")
    cols=cols[1:length(folders)]
    
    # go throught each folder for plotting log and linear axes
    par(mfrow = c(1,1))
    
    # Linear Scale---------------------------
    png(paste0(outputfolder,"KM_plots_linear.png"),width=1200,height=800)
    par(mar=c(5.1,5.1,4.1,2.1))   
    par(oma=c(0,0,0,0))  
    par(cex.lab=2)  	     
    par(cex.main=2)        
    par(cex.axis=2)
    
    for (i in 1:length(folders)){
      if ( i==1 ){
        actual = all_tab[[i]]
        actual$SurvObj <- with(actual, Surv(TIME, DV %in% c(1,3)))
        km.as.one <- survfit(SurvObj ~ 1, data = actual )
        plot(km.as.one, mark.time=TRUE, mark=1, ylim = c(0,1),
             ... )
      }
      lines( all_tab[[i+1]], col=cols[i], lty = i+1,lwd=2 )
      if (i==length(folders)){
        legend("bottomleft", model_names, col = cols,
               lty = 2:(length(model_names)+1),
               lwd=2,pch = NA,
               merge = TRUE, bg = NA,
               cex=2
        )
      }
    }
    dev.off()
    
    # log scale---------------------------
    png(paste0(outputfolder,"KM_plots_log.png"),width=1200,height=800)
    par(mar=c(5.1,5.1,4.1,2.1))   
    par(oma=c(0,0,0,0))  
    par(cex.lab=2)  	     
    par(cex.main=2)        
    par(cex.axis=2)
    
    for (i in 1:length(folders)){
      if ( i==1 ){
        actual = all_tab[[i]]
        actual$SurvObj <- with(actual, Surv(TIME, DV %in% c(1,3)))
        km.as.one <- survfit(SurvObj ~ 1, data = actual )
        plot(km.as.one, mark.time=TRUE, mark=1, ylim = c(0,1),
             log = "x",
             ... )
      }
      lines( all_tab[[i+1]], col=cols[i], lty = i+1,lwd=2 )
      
      if (i==length(folders)){
        legend("bottomleft", model_names, col = cols,
               lty = 2:(length(model_names)+1),
               lwd=2,pch = NA,
               merge = TRUE, bg = NA,
               cex=2
        )
      }
      
    }
    dev.off()
    #########################################
    #########################################
  }
}
##############################################################
###### END OF FUNCTION Models_summary ########################
##############################################################


##############################################################
### FUNCTION FOR CALCULATING HAZARD FUNCTION #################
##############################################################
### LIST OF FUNCTION ARGUMENTS
#
# t_vals                    # time values at which to output value of the hazard function
# parameters:               # vector of parameter values which define the hazard function
#
# Required entries within 'parameters':
# HazardFunction            # text defining hazard function type. Six possible values: exponential / gompertz / weibull / lomax / log_logistic / log_normal
# lambda                    # lambda value for all six possible hazard function types.
# alpha                     # alpha value for gompertz / weibull / lomax / log_logistic / log_normal hazard functions. Set to NA if HazardFunction='exponential'.
# CovMod                    # a value of 1 or 0 which defines if the hazard function is influenced by covariate(s). 0 means no covariate influence, 1 means there is a covariate influence.
#
# Definitions of the 6 possible hazard functions.
# HazardFunction='exponential':            h(t) = lambda
# HazardFunction='gompertz':               h(t) = lambda*exp(alpha*t)
# HazardFunction='weibull':                h(t) = lambda*exp(alpha*ln(t))
# HazardFunction='lomax':                  h(t) = lambda/(alpha+t)
# HazardFunction='log_logistic':           h(t) = (alpha*(lambda^alpha)*t^(alpha-1)/(1+(lambda^alpha)*t^alpha))
# HazardFunction='log_normal':             h(t) = pdf/(1-phi(log(t/alpha)/lambda))
#                                                 Where pdf is the probability density function of the log-normal distribution,
#                                                 parameterised as Lognormal(log(alpha),lambda) and phi is the 
#                                                 cumulative distribution function of the normal distribution.
# HazardFunction='log_cauchy':             h(t) = pdf/(1-cdf)
#                                                 Where pdf = 1/(t*pi)*(lambda/(log(t/alpha)^2 + lambda^2))
#                                                 and cdf = (1/pi)*arctan(log(t_vals/alpha)/lambda) + 0.5
#
# NOTE: In order to simulate with covariate influence(s) on the hazard function, the function cov_func() needs to have been defined. This is done by the sim_events() function.
##############################################################
Hazard=function(t_vals,parameters){
  # extract parameter values
  HazFunc=parameters['HazardFunction']
  CovInfluence=as.numeric(parameters['CovMod'])   
  if(HazFunc=='exponential'){
    Lambda=as.numeric(parameters['lambda'])
  }
  if(HazFunc%in%c('gompertz','weibull','lomax','log_logistic','log_normal','log_cauchy')){
    Lambda=as.numeric(parameters['lambda'])
    Alpha=as.numeric(parameters['alpha'] )  
  }
  
  # build the time-dependent haz_fac if a covariate model is required
  if(CovInfluence==0){
    haz_fac=1
  }else{ 
    # calculate value of time-varying covariate for current value of t (using the cov_func() function, as calculated within the sim_events() function)
     if(is.null(cov_func)==TRUE){
       stop('Error in function Hazard(): Covariate model was requested but cov_func() does not exist')
     }else{
       haz_fac=cov_func(t_vals)
     }
  }

  # calculate the hazard function   
  
  DEL = 1E-08              # Some of the models require a small number to be added to time for stability
  
  if(HazFunc=='exponential'){
    Res=haz_fac*Lambda
  }
  if(HazFunc=='gompertz'){
    Res=haz_fac*Lambda*exp(Alpha*t_vals)
  }
  if(HazFunc=='weibull'){
    Res=haz_fac*Lambda*exp((Alpha-1)*log(t_vals+DEL))
  }
  if(HazFunc=='lomax'){
    Res=haz_fac*Lambda/(Alpha+t_vals)
  }
  if(HazFunc=='log_logistic'){
    Res=haz_fac*(Alpha*(Lambda^Alpha)*(t_vals+DEL)^(Alpha-1)/(1+(Lambda^Alpha)*(t_vals+DEL)^Alpha))
  }
  if(HazFunc=='log_normal'){
    t_vals=t_vals+DEL
    num = log(t_vals/Alpha)
    pdfunc = 1/(Lambda*sqrt(2*pi)*t_vals)*exp(-num^2/(2*Lambda^2))
    Res=haz_fac*pdfunc/(1-pnorm(num/Lambda))
  }
  if(HazFunc=='log_cauchy'){
    t_vals=t_vals+DEL
    num = log(t_vals/Alpha)/Lambda
    pdfunc = 1/(t_vals*pi)*(Lambda/(log(t_vals/Alpha)^2 + Lambda^2))
    cdfunc = (1/pi)*atan(num) + 0.5
    Res=haz_fac*pdfunc/(1-cdfunc)
  }
  return(Res)
}
##############################################################
###### END OF FUNCTION Hazard ################################
##############################################################


###################################################################
### FUNCTION FOR INTEGRATING HAZARD TO GIVE CUMULATIVE HAZARD #####
###################################################################
### LIST OF FUNCTION ARGUMENTS
#
#  t                        # time t as sent by ode function of deSolve library
#  y                        # initial value for compartment A1 [should always be zero]      
#  parms                    # vector of parameter values that will be passed on to the Hazard function defined above [see above for definition of content of parms]
###################################################################
CumHaz=function(t,y,parms){
  with(as.list(c(y,parms)), {
    dA1 = Hazard(t,parms)
    return(list(dA1))
  })
}
##################################################################
###### END OF FUNCTION CumHaz ####################################
##################################################################


###################################################################
### WRAPPER FUNCTION FOR INTEGRATING HAZARD TO GIVE CUMULATIVE ####
### HAZARD AND SURVIVAL FUNCTION ##################################
###################################################################
### LIST OF FUNCTION ARGUMENTS
#
#  parms                    # vector of parameter values that will be passed on to the Hazard function defined above [see above for definition of content of parms]
#  maxtime                  # maximum time for integration
#  timeint                  # time interval for result output
#
# NOTE: this function is not for use with hazard functions that include covariate influences. If parms contains CovMod of 1, then CovMod will be set to 0.
###################################################################
IntHaz=function(parms=NA,maxtime=1000,timeint=1){
  if(as.character(parms['CovMod'])=='1'){
     parms['CovMod']='0'
  }
  tvals=seq(0,ceiling(maxtime),by=timeint)
  res=suppressWarnings(ode(y=c(A1=0),times=tvals,func=CumHaz,parms=parms))
  res=data.frame(res)
  names(res)[2]='CHAZ'
  res$SURV=exp(-res$CHAZ)
  return(res)
}
##################################################################
###### END OF FUNCTION IntHaz ####################################
##################################################################



##############################################################
### FUNCTION FOR PLOTTING KAPLAN-MEIER PLOTS #################
##############################################################
### LIST OF FUNCTION ARGUMENTS
#
#  event_data_frame         # data frame containing event times and event types (1=event, 0=censoring)
#  time_col                 # column name for event times       
#  event_col                # columbn name for event types
#  conf_interval_pct        # % value for KM confidence interval. Set to FALSE if KM confidence intervals are not required
#  x_range                  # display range of x axis eg c(0,1000)
#  y_range                  # display range of y axis eg c(0,1)
#  x_label                  # label for x axis
#  y_label                  # label for y axis
#  line_width               # line width for the KM curve and confidence bands (typically 1 or 2)
#  KM_colour                # colour for the KM curve. Set this to NA if you do not want the KM curve added to the plot
#  conf_int_colour          # colour for the dashed line KM confidence intervals
#  mark_censored_events     # set to TRUE/FALSE
#  add_to_existing_plot     # normally set to FALSE. Set to TRUE for adding KM curve on top of existing plot of simulation confidence bands so the KM curve can be seen better
##############################################################
plot_KM=function(event_data_frame,time_col='TIME',event_col='DV',conf_interval_pct=95,x_range=c(0,1000),y_range=c(0,1),x_label='Time (days)',
                 y_label='Probability of event at time > t',line_width=2,KM_colour='black',conf_int_colour='black',
                 mark_censored_events=TRUE,add_to_existing_plot=FALSE,...){ 
  
  s1=Surv(event_data_frame[,time_col],event_data_frame[,event_col])
  if(conf_interval_pct!=FALSE){
    CI=conf_interval_pct/100
    conf_int=CI
    use_conf_int=TRUE
  } else {
    use_conf_int=FALSE
    conf_int=0.95
    conf_int_colour=KM_colour
  } 
  km1=survfit(s1~1,conf.int=conf_int)
  if(add_to_existing_plot==TRUE){
    par(new=TRUE)
    ax=FALSE
    x_label=''
    y_label=''
  } else {
    ax=TRUE
  }
  co=c(KM_colour,conf_int_colour,conf_int_colour)
  plot(km1,mark.time=mark_censored_events,mark=0,conf.int=use_conf_int,xlim=x_range,ylim=y_range,xlab=x_label,ylab=y_label,lwd=line_width,col=co,axes=ax,...)
}
##############################################################
###### END OF FUNCTION plot_KM ###############################
##############################################################


#################################################################################################################
### function to simulate individual survival and hazard curves and an event for each patient in a study #########
#################################################################################################################
###
### Function Inputs:
### Patient_Info:   dataframe with 2 required columns: 1) ID column contains patient IDs. 2) a column that indicates the time length for which the hazard function should be integrated for each patient
### Int_Time_Col:   name of the column within Patient_Info that contains the time length for integration
### Parameters:     a vector of named elements which define parameters of the hazard function
### ls_times:       for defining covariate influences: a list object with time data for each individual.
### ls_cov:        for defining covariate influences: a list object which quantifies covariate influence for each individual.
### Method:         integration method to be used by desolve, default lsoda, if stiff, try lsode or bdf
###
### NOTE: THE HAZARD FUNCTION WILL ONLY BE INTEGRATED UNTIL THE VALUE WITHIN Int_Time_Col FOR EACH PATIENT. 
###       IF NO EVENT SIMULATED BY THIS TIME THEN EVENT CENSORED AT THAT TIME.
###       IF VALUE IN Int_Time_Col is > THE MAXIMUM TIME WHERE COVARIATE DATA IS AVAILABLE THEN THE FINAL VALUE
###       OF THE COVARIATE IS CARRIED FORWARDS IN THE INTEGRATION OF THE HAZARD FUNCTION
###
### Function Outputs:
### A list object containing:
###    $sim_events. A data frame containing a single set of simulated event times, with these columns:
### "event_time" contains the time of the simulated event
### "event" contains 1 for events and 0 for censoring
### "patient_ID" contains patient identifiers
### "integration_time" contains the end time of the hazard function integration
###
###    $surv_haz. A list object containing data frames of individual survival and hazard curves for each patient. 
### Each data frame within $surv_haz is named by the patient identifier. Each data frame contains these columns:
### "time" contains the time
### "CHAZ" contains individual cumulative hazard
### "Survival" contains individual survival
### "Hazard" contains individual hazard
### "ID" contains patient identifiers
### "integration_time" contains the end time of the hazard function integration
##############################################################
sim_events=function(Patient_Info=NA,Int_Time_Col=NA,Parameters=NA,ls_times=NA,ls_cov=NA,Method="lsoda"){
   
  patient_IDs=as.character(Patient_Info$ID)
  
  num_patients=nrow(Patient_Info)
  
  model=Parameters['HazardFunction']
  CovInfluence=as.numeric(Parameters['CovMod'])   
  
  event=rep(0,num_patients)
  event_time=rep(-9,num_patients)
  patient_number=rep(-9,num_patients)
  U_vals=runif(n=num_patients,min=0,max=1)
  surv_haz_curves=vector(mode='list',length=num_patients)
  names(surv_haz_curves)=patient_IDs
  Y_ini=c(A1=0)
  par(mfrow=c(1,1))
  for (i in 1:num_patients){
    
    print(paste0('Calculating individual survival curve: ',i,'/',num_patients,':  Patient ID=',Patient_Info[i,1]))
    
    # assign random number in range 0 to 1, to current patient
    U_val=U_vals[i]
    
    # set up cov_func() for use by Hazard() when hazard function is integrated for current individual
    # make sure any existing instance of cov_func is removed
    cov_func=NULL
    cov_func<<-NULL
    if(CovInfluence==1){
      # build data frame of time-dependent covariate data for current patient
      SID=patient_IDs[i]
      cov_data=data.frame(time=ls_times[[SID]],Cov_Values=ls_cov[[SID]])
      # generate a function for interpolation of covariate data at any value of time for the current patient
      cov_func<<-approxfun(cov_data,method='constant',rule=2,f=1)
    }

    # now integrate the hazard function (up to time given by value in Int_Time_Col)
    T_vals=seq(from=0,to=Patient_Info[i,Int_Time_Col],by=1)
    if(Patient_Info[i,Int_Time_Col]%%1>0){
      dec_val=Patient_Info[i,Int_Time_Col]%%1
      T_vals=c(T_vals,max(T_vals)+dec_val)
    }
    
    surv=ode(y=Y_ini,times=T_vals,func=CumHaz,parms=Parameters,method=Method)
    surv=data.frame(surv)
    names(surv)[2]='CHAZ'
    surv_res=surv
    surv=exp(-surv$CHAZ)
    surv_res$Survival=surv
    
    # save the hazard and survival curves to list objects
    xz=Hazard(T_vals,Parameters)
    xz=data.frame(xz)
    names(xz)='Hazard'
    surv_res$Hazard=xz$Hazard
    surv_res$ID=Patient_Info[i,1]
    surv_res$integration_time=Patient_Info[i,Int_Time_Col]
    surv_haz_curves[[i]]=surv_res
    
    # find point in time where survival curve drops below current random number U_val
    # if no event is simulated within time=0 and time=Int_Time_Col then patient is censored at Int_Time_Col
    right_censor_time=max(T_vals)
    q=surv<=U_val
    test=T_vals[q]
    sim_event_time=right_censor_time
    event_type=0
    if(length(test)>0){
      sim_event_time=test[1]
      event_type=1
    } 
    event_time[i]=sim_event_time
    event[i]=event_type
    patient_number[i]=patient_IDs[i]
  }
  patient_ID=patient_number
  sim_events=data.frame(event_time,event,patient_ID)
  sim_events$integration_time=Patient_Info[,Int_Time_Col]
  res=list(sim_events=sim_events,surv_haz=surv_haz_curves)
  
  print(paste0('Individual survival curves calculated using ',model,' hazard function'))
  return(res)
}
##############################################################
####### END OF FUNCTION sim_events ###########################
##############################################################


########################################################################################
### function to produce replicate simulations of events for each patient in a study
########################################################################################
### Function Inputs:
### Patient_Info:   dataframe with 1 required columns: ID column contains patient IDs.
### num_sim_repeats:  number of times the entire simulation should be repeated
### survival_curves:  list of individual survival curves for events from output of sim_events function
###
### Function Outputs:
### A list object containing:
###    $events. A list object containing one data frame for each simulation repeat.
###    Each data frame has these columns: 
###                                       patient_ID
###                                       event_time
###                                       event (1=event, 0=right censored)
###                                       max_cov_time, which is maximum time at which covariate data is available. If no event is
###                                       simulated up until max_cov time then right censoring is at max_cov_time
##########################################################################################
repeated_sim_events=function(Patient_Info=NA,num_sim_repeats=NA,survival_curves=NA){
  
  num_patients=nrow(Patient_Info)
  
  # first define a list which will contain a data frame of event times and types for each simulation repeat
  rep_event_list=vector(mode='list',length=num_sim_repeats)
  rep_dropout_list=vector(mode='list',length=num_sim_repeats)
  
  # generate events
  for(i in 1:num_sim_repeats){
    print(paste('Replicate: ',i))
    U_vals=runif(n=num_patients,min=0,max=1)
    event=rep(0,num_patients)
    event_time=rep(-9,num_patients)
    patient_ID=rep(-9,num_patients)
    integration_time=rep(-9,num_patients)
    for(z in 1:num_patients){
      # find point in time where survival curve drops below current random number U_val
      T_vals=survival_curves[[z]]['time']
      right_censor_time=max(T_vals)
      surv=survival_curves[[z]]['Survival']
      q=surv<=U_vals[z]
      test=T_vals[q]
      sim_event_time=right_censor_time
      event_type=0
      if(length(test)>0){
        sim_event_time=test[1]
        event_type=1
      } 
      event_time[z]=sim_event_time
      event[z]=event_type
      patient_ID[z]=Patient_Info[z,'ID']
      integration_time[z]=survival_curves[[z]][,'integration_time'][1]
    }
    rep_event_list[[i]]=data.frame(patient_ID,event_time,event,integration_time)
  }
  
  res=list(events=rep_event_list)
  return(res)
  print('Finished generating repeat sets of events.')
}
##############################################################
####### END OF FUNCTION repeated_sim_events ##################
##############################################################


###########################################################################################
### FUNCTION TO ADJUST SIMULATED EVENT TIMES ACCORDING TO ANALYSIS OF COMPETING EVENTS  ###
###########################################################################################
### Function Inputs:
### event_name:        a name for the event type simulated within "events" list object
### events:            a list object (output of repeated_sim_events()) containing many data frames of simulated events. 
### competing_events:  a list of list objects, each sub-list being the output of repeated_sim_events() and containing simulated events that compete with simulated events within "events".
###
### Function Outputs:
### A list object containing:
###    $events. A list object containing one data frame for each simulation repeat.
###    Each data frame includes these columns (among others): 
###                                       patient_ID
###                                       event_time          original event time before competing event analysis
###                                       event               original event type: 1=event of type event_name, 0=right censored from event_name    
###                                       new_event_time      event time after ajustment from competing event analysis
###                                       new_event           new event type: 1=event of type event_name, 0=right censored from event_name 
###                                       new_event_code      gives a text description of the simulated event, including a reason for censoring
###
### NOTE: If any of the the simulations contain "adjusted_event_time" (from function adjust_events_to_appointment_times()) then that will be used in place of "event_time".
##########################################################################################
CompeteEvents=function(event_name=NA,events=NA,competing_events=NA){
  comp_names=names(competing_events)   # get names of the competing events
  sim_n=length(events$events)          # get number of simulation replicates
  # check that all event objects have same number of simulation replicates
  for(i in 1:length(competing_events)){
    if(length(competing_events[[comp_names[i]]]$events)!=sim_n){
      msg=paste0('ERROR: Not all event objects have same number of simulation replicates')
      print(msg)
      stop()
    }
  }
  
  # if any of the simulation objects contain "adjusted_event_time" then we need to use that,
  # so copy "adjusted_event_time" to "event_time"
  if('adjusted_event_time'%in%names(events$events[[1]])==TRUE){
    for(i in 1:sim_n){
      events$events[[i]]$event_time=events$events[[i]]$adjusted_event_time
    }
  }
  for(i in 1:length(competing_events)){
    if('adjusted_event_time'%in%names(competing_events[[comp_names[i]]]$events[[1]])==TRUE){
      for(j in 1:sim_n){
        competing_events[[comp_names[i]]]$events[[j]]$event_time=competing_events[[comp_names[i]]]$events[[j]]$adjusted_event_time
      }
    }
  }
  
  # set up names for the competing event columns that will be copied to "events" data frames
  comp_ev_cols=character(length(comp_names))
  comp_ev_t_cols=character(length(comp_names))
  ev_cols=character(length(comp_names))
  for(i in 1:length(comp_names)){
    comp_ev_cols[i]=paste0('event_',comp_names[i])
    comp_ev_t_cols[i]=paste0('time_',comp_names[i])
  }
  all_t_cols=c('event_time',comp_ev_t_cols)
  all_ev_cols=c(event_name,comp_names)
  
  # now loop through each simulation replicate and do the competing event analysis
  for(i in 1:sim_n){
    # copy the competing events and event times to 'events'
    for(k in 1:length(comp_names)){
      events$events[[i]]$new_event_time=-5
      events$events[[i]]$new_event=-5
      events$events[[i]]$new_event_code=-5
      events$events[[i]][,comp_ev_t_cols[k]]=competing_events[[k]]$events[[i]]$event_time
      events$events[[i]][,comp_ev_cols[k]]=competing_events[[k]]$events[[i]]$event
    }
    
    # compare time of the various simulated events
    z=events$events[[i]]
    q=events$events[[i]][,c(all_t_cols)]
    q[,]=NA
    q[z$event==1,'event_time']=z[z$event==1,'event_time']
    for(x in 1:length(comp_names)){
      q[z[,comp_ev_cols[x]]==1,comp_ev_t_cols[x]]=z[z[,comp_ev_cols[x]]==1,comp_ev_t_cols[x]]
    }
    q1=q
    q1[is.na(q1)]=10*max(q1,na.rm=TRUE)
    q$min_time=suppressWarnings(apply(q,1,min,na.rm=TRUE))
    q$res=apply(q1,1,which.min)
    q$res[q$min_time==Inf]=NA
    q$res1=all_ev_cols[q$res]
    events$events[[i]]$new_event_time=q$min_time
    events$events[[i]]$new_event_code=q$res1
    
    # we are left with some events that have not been triggered before integration_time
    events$events[[i]]$new_event_time[is.na(events$events[[i]]$new_event_code)==TRUE]=events$events[[i]]$integration_time[is.na(events$events[[i]]$new_event_code)==TRUE]
    events$events[[i]]$new_event_code[is.na(events$events[[i]]$new_event_code)==TRUE]='no_sim_events_before_integration_end_time'
    # finally set the new_event
    events$events[[i]]$new_event=0
    events$events[[i]]$new_event[events$events[[i]]$new_event_code==event_name]=1
    
  }
  return(events)
}
##############################################################
####### END OF FUNCTION CompeteEvents ########################
##############################################################


########################################################################################
### FUCNTION TO CARRY OUT CALCULATION OF EVENT VPC #####################################
########################################################################################
### Function Inputs:
### Patient_Info:           dataframe with 2 required columns; ID column contains patient IDs, and covariate column used for stratification as defined by the 'stratify' argument
### obs_events:             dataframe containing the observed events
### obs_event_cols:         vector containing column header for time and event (0/1) within obs_events. Order has to be time then event.
### sim_events:             list object containing replicate simulated data. Object is the output of either repeated_sim_events() or CompeteEvents().
### compete:                TRUE/FALSE. If sim_events is the output from repeated_sim_events() then set compete=FALSE. If sim_events is the output from CompeteEvents() then set compete=TRUE or FALSE to investigate the effect of analysing competing events.
### bins:                   time points for counting of cumulative events. Events go in a bin if event_time<=bin value.
### conf_int:               confidence interval width. Defaults to 90%.
### stratify:               name of a stratification column within Patient_Info. Set to FALSE if not used.
###
### Function Outputs:
### A list object containing:
###  $values:               dataframe containing VPC results as cumulative counts up to (and including) each time within 'bins'
###  $values_noncumulative: dataframe containing VPC results as discrete counts between each time within 'bins' (strictly: greater than bins[n] and less than or equal to bins[n+1])
###  $stratification:       value of "stratify" when function was exectuted
###  $num_strat_levels:     number of levels within the stratification covariate
###  $strat_levels:         values of the stratification covariate
###  $strat_counts:         counts of patients within each level of the stratification covariate
###  $confidence_interval:  confidence interval used during calculation 
###
### NOTE: If the simulations contain "adjusted_event_time" (from function adjust_events_to_appointment_times()) then that will be used in place of "event_time".
##########################################################################################
calcEventVPC=function(Patient_Info=NA,obs_events=NA,obs_event_cols=c('TIME','DV'),sim_events=NA,
                      compete,bins=c(100,200,300,400,500,600,700,800,900,1000),conf_int=90,stratify=FALSE){
  
  n_rep=length(sim_events$events)
  CI_lower=(100-conf_int)/2
  CI_upper=100-(100-conf_int)/2
  
  if(compete==TRUE){
    sim_event_cols=c('new_event_time','new_event')
  }else{
    sim_event_cols=c('event_time','event')
  }
  
  print('calculating VPC...')
  
  # if the simulation object contains "adjusted_event_time" then we need to use that,
  # so copy "adjusted_event_time" to "event_time"
  if('adjusted_event_time'%in%names(sim_events$events[[1]])==TRUE){
    for(i in 1:n_rep){
      sim_events$events[[i]]$event_time=sim_events$events[[i]]$adjusted_event_time
    }
  }
  
  # set up the data frame that will contain calculation results
  strat_levels=0
  n_strat_levels=0
  st_counts=nrow(Patient_Info)
  if(stratify!=FALSE){
    bb=data.frame(table(Patient_Info[,stratify]))
    bb=bb[order(bb$Var1),]
    st_counts=bb$Freq
    strat_levels=as.character(bb$Var1)
    n_strat_levels=length(strat_levels)
  }
  res=data.frame(matrix(data=NA,nrow=length(bins),ncol=4*n_strat_levels+5))
  ns=c('obs_events','median_sim_events',paste0('LL_CI_sim_events'),paste0('UL_CI_sim_events'))
  if(n_strat_levels>1){
    cnms=c('time',ns)
    for(i in 1:n_strat_levels){
      stb=paste0(strat_levels[i],'_')
      ns1=paste0(stb,ns)
      cnms=c(cnms,ns1)
    }
  }else{
    cnms=c('time',ns)
  }
  names(res)=cnms
  res$time=bins
  res1=res
  
  # analyze observed events (unstratified)
  bins1=c(0,bins)
  for(i in 1:length(bins)){
    d=obs_events[,obs_event_cols[2]][obs_events[,obs_event_cols[1]]<=bins[i]]
    res$obs_events[i]=sum(d)
    if(i==1){
       d=obs_events[,obs_event_cols[2]][obs_events[,obs_event_cols[1]]<=bins1[i+1] & obs_events[,obs_event_cols[1]]>=bins1[i]]
       res1$obs_events[i]=sum(d)
    }else{
       d=obs_events[,obs_event_cols[2]][obs_events[,obs_event_cols[1]]<=bins1[i+1] & obs_events[,obs_event_cols[1]]>bins1[i]]
       res1$obs_events[i]=sum(d)
    }
  }
  # analyze observed events (stratified)
  if(stratify!=FALSE){
    for(i in 1:n_strat_levels){
       dat=obs_events[obs_events[,stratify]==strat_levels[i],]
       for(k in 1:length(bins)){
         d=dat[,obs_event_cols[2]][dat[,obs_event_cols[1]]<=bins[k]]
         chead=paste0(strat_levels[i],'_obs_events')
         res[k,chead]=sum(d)
         if(k==1){
            d=dat[,obs_event_cols[2]][dat[,obs_event_cols[1]]<=bins1[k+1] & dat[,obs_event_cols[1]]>=bins1[k] ]
            chead=paste0(strat_levels[i],'_obs_events')
            res1[k,chead]=sum(d)
         }else{
            d=dat[,obs_event_cols[2]][dat[,obs_event_cols[1]]<=bins1[k+1] & dat[,obs_event_cols[1]]>bins1[k] ]
            chead=paste0(strat_levels[i],'_obs_events')
            res1[k,chead]=sum(d)
         }
       }
    } 
  }
  
  # analyze simulated events (unstratified)
  vals=data.frame(matrix(data=NA,nrow=length(bins),ncol=n_rep))
  vals1=vals
  for(k in 1:n_rep){
    dat=sim_events$events[[k]]
      for(i in 1:length(bins)){
        d=dat[,sim_event_cols[2]][dat[,sim_event_cols[1]]<=bins[i]]
        vals[i,k]=sum(d)
        if(i==1){
           d=dat[,sim_event_cols[2]][dat[,sim_event_cols[1]]<=bins1[i+1] & dat[,sim_event_cols[1]]>=bins1[i]]
           vals1[i,k]=sum(d)
        }else{
           d=dat[,sim_event_cols[2]][dat[,sim_event_cols[1]]<=bins1[i+1] & dat[,sim_event_cols[1]]>bins1[i]]
           vals1[i,k]=sum(d)
        }
      }
  }
  # calculate result entries from data within vals / vals1
  for(i in 1:length(bins)){
     res$median_sim_events[i]=median(as.numeric(vals[i,]))
     title=paste0('LL_CI_sim_events')
     res[i,title]=quantile(as.numeric(vals[i,]),probs=CI_lower/100)
     title=paste0('UL_CI_sim_events')
     res[i,title]=quantile(as.numeric(vals[i,]),probs=CI_upper/100)
     res1$median_sim_events[i]=median(as.numeric(vals1[i,]))
     title=paste0('LL_CI_sim_events')
     res1[i,title]=quantile(as.numeric(vals1[i,]),probs=CI_lower/100)
     title=paste0('UL_CI_sim_events')
     res1[i,title]=quantile(as.numeric(vals1[i,]),probs=CI_upper/100)
  }
  
  # analyze observed events (stratified)
  if(stratify!=FALSE){
    for(z in 1:n_strat_levels){
      vals=data.frame(matrix(data=NA,nrow=length(bins),ncol=n_rep))
      for(k in 1:n_rep){
        dat=sim_events$events[[k]]
        dat[,stratify]=Patient_Info[,stratify]    # copy the stratification covariate from Patient Info
        dat=dat[dat[,stratify]==strat_levels[z],]
        for(i in 1:length(bins)){
          d=dat[,sim_event_cols[2]][dat[,sim_event_cols[1]]<=bins[i]]
          vals[i,k]=sum(d)
          if(i==1){
             d=dat[,sim_event_cols[2]][dat[,sim_event_cols[1]]<=bins1[i+1] & dat[,sim_event_cols[1]]>=bins1[i]]
             vals1[i,k]=sum(d)
          }else{
             d=dat[,sim_event_cols[2]][dat[,sim_event_cols[1]]<=bins1[i+1] & dat[,sim_event_cols[1]]>bins1[i]]
             vals1[i,k]=sum(d)
          }
        }
      }
      # calculate result entries from data within vals
      for(h in 1:length(bins)){
        title=paste0(strat_levels[z],'_median_sim_events')
        res[h,title]=median(as.numeric(vals[h,]))
        title=paste0(strat_levels[z],'_LL_CI_sim_events')
        res[h,title]=quantile(as.numeric(vals[h,]),probs=CI_lower/100)
        title=paste0(strat_levels[z],'_UL_CI_sim_events')
        res[h,title]=quantile(as.numeric(vals[h,]),probs=CI_upper/100)
        title=paste0(strat_levels[z],'_median_sim_events')
        res1[h,title]=median(as.numeric(vals1[h,]))
        title=paste0(strat_levels[z],'_LL_CI_sim_events')
        res1[h,title]=quantile(as.numeric(vals1[h,]),probs=CI_lower/100)
        title=paste0(strat_levels[z],'_UL_CI_sim_events')
        res1[h,title]=quantile(as.numeric(vals1[h,]),probs=CI_upper/100)
      }
    }
  }
  print('finished calculating VPC')
  if(stratify!=FALSE){
    stratif=stratify
  }else{
    stratif='none'
  }
  
  result=list(values=res,values_noncumulative=res1,stratification=stratif,num_strat_levels=n_strat_levels,strat_values=strat_levels,strat_counts=st_counts,confidence_interval=conf_int)
  return(result)
}
##############################################################
####### END OF FUNCTION calcEventVPC #########################
##############################################################


########################################################################################
### FUNCTION TO PLOT EVENT VPC #########################################################
########################################################################################
### Function Inputs:
### VPCdata:                VPC calculation results from output of calcEventVPC()
### stratify:               TRUE/FALSE. Instruction to produce a stratified VPC
### type:                   'cumulative'=display cumulative events on y-axis. 'discrete'=display discrete counts on y-axis.
### display_percent:        TRUE/FALSE. Set to TRUE to display cumulative events count as percentage of patient count
### x_label:                x-axis label for VPC plot, defaults to 'Time (days)'
### y_label:                y-axis label for VPC plot, defaults to 'Cumulative event count'
### bin_width:              width of each display bin, in time units
### strat_spacing:          spacing between each level of stratification, in units of time
### col_scheme:             colour scheme, a vector of three colours, default=c('red','black','gray90','gray80','gray70','gray60','gray50'). 1st value=observed events, 2nd=median of sim events, 3rd onwards=conf interval shading
### add_strat_lines:        TRUE/FALSE. Instruction to add dotted vertical lines between stratification results
### legend_position:        legend position direction: "bottomright", "bottom","bottomleft", "left", "topleft", "top", "topright","right", "center"
### legend_size:            legend size. Default=1.2
###
##########################################################################################
plotEventVPC=function(VPCdata=NA,stratify=FALSE,type='cumulative',display_percent=FALSE,x_label='Time (days)',y_label='Cumulative event count',
                      bin_width=15,strat_spacing=5,col_scheme=c('red','black','gray90','gray80','gray70','gray60','gray50'),add_strat_lines=TRUE,
                      legend_position='bottomright',legend_size=1,...){
  
  xrange=c(0,max(VPCdata$values$time))   # set x-axis display range
  
  if(legend_position!=FALSE){
    if (!legend_position %in% c("bottomright", "bottom","bottomleft", "left", "topleft", "top", "topright","right", "center")) {
      stop("Legend position not recognised, pleasse use one of the following: bottomright, bottom, bottomleft, left, topleft, top, topright, right, center")
    }
  }
  
  fillcol=col_scheme[3]
  sim_med=col_scheme[2]
  obs_lines=col_scheme[1]
  
  if(stratify==TRUE){
    if(VPCdata$stratification=='none'){
      stop('Cannot plot stratified VPC because stratification calculations have not been performed.')
    }
  }
  
  if(display_percent==TRUE & y_label=='Cumulative event count'){
    y_label='% patients with event at time < t'
  }
  
  # set up stratification spacing
  if(stratify==TRUE){
    display_width=VPCdata$num_strat_levels*bin_width+(VPCdata$num_strat_levels-1)*strat_spacing
  }else{
    display_width=bin_width
  }
  start_offset=-display_width/(VPCdata$num_strat_levels*2)  # start drawing line from time=time value - start_offset
  
  par(mfrow=c(1,1)) 
  par(mar=c(4.5,4.5,2,2))            
  par(oma=c(0,0,0,0))
  par(cex.lab=1.5)  	
  par(cex.main=1.5)        
  par(cex.axis=1.5)
  
  
  if(type%in%c('cumulative')){
     d=VPCdata$values
     
     # convert to percentages if necessary
     if(display_percent==TRUE){
       if(stratify==FALSE){
          cnt=sum(VPCdata$strat_counts)
          d$obs_events=100*d$obs_events/cnt
          d$median_sim_events=100*d$median_sim_events/cnt
          d$LL_CI_sim_events=100*d$LL_CI_sim_events/cnt
          d$UL_CI_sim_events=100*d$UL_CI_sim_events/cnt
       }else{
         for(i in 1:VPCdata$num_strat_levels){
           s_nm=paste0(VPCdata$strat_values[i],'_')
           sim_med_name=paste0(s_nm,'median_sim_events')
           sim_LL_CI_name=paste0(s_nm,'LL_CI_sim_events')
           sim_UL_CI_name=paste0(s_nm,'UL_CI_sim_events')
           obs_name=paste0(s_nm,'obs_events')
           d[,obs_name]=100*d[,obs_name]/VPCdata$strat_counts[i]
           d[,sim_med_name]=100*d[,sim_med_name]/VPCdata$strat_counts[i]        
           d[,sim_LL_CI_name]=100*d[,sim_LL_CI_name]/VPCdata$strat_counts[i]
           d[,sim_UL_CI_name]=100*d[,sim_UL_CI_name]/VPCdata$strat_counts[i]
         }
       }
     }
     
     n_bins=nrow(d)
     if(stratify==TRUE){
       d$obs_events=NULL
       d$median_sim_events=NULL
       d$LL_CI_sim_events=NULL
       d$UL_CI_sim_events=NULL
       plot(x=xrange,y=c(0,max(d[nrow(d),-1])),type='n',xlab=x_label,ylab=y_label,...)
     }else{
       plot(x=xrange,y=c(0,max(d[nrow(d),c('obs_events','median_sim_events','LL_CI_sim_events','UL_CI_sim_events')])),type='n',xlab=x_label,ylab=y_label,...)
     }
     
     if(stratify==FALSE){
        for(i in 1:n_bins){
          polygon(x=c(d$time[i]-bin_width/2,d$time[i]+bin_width/2,d$time[i]+bin_width/2,d$time[i]-bin_width/2),y=c(d$LL_CI_sim_events[i],d$LL_CI_sim_events[i],d$UL_CI_sim_events[i],d$UL_CI_sim_events[i]),col=fillcol,border=fillcol)
          lines(x=c(d$time[i]-bin_width/2,d$time[i]+bin_width/2),y=c(d$obs_events[i],d$obs_events[i]),lwd=2,col=obs_lines) 
          lines(x=c(d$time[i]-bin_width/2,d$time[i]+bin_width/2),y=c(d$median_sim_events[i],d$median_sim_events[i]),lwd=1,col=sim_med) 
        }
     }
     
     if(stratify==TRUE){
       for(k in 1:VPCdata$num_strat_levels){
         s_level=VPCdata$strat_values[k]
         s_name=paste0(s_level,'_')
         
         if(k>1){
           start_offset=start_offset+bin_width+strat_spacing
         }
         add_vert=FALSE
         if(k<VPCdata$num_strat_levels & add_strat_lines==TRUE){
           add_vert=TRUE
           next_s_level=VPCdata$strat_values[k+1]
           next_s_name=paste0(next_s_level,'_')
         }
           for(i in 1:n_bins){
             ctime=d$time[i]
             x_start=ctime+start_offset
             sim_med_name=paste0(s_name,'median_sim_events')
             sim_LL_CI_name=paste0(s_name,'LL_CI_sim_events')
             sim_UL_CI_name=paste0(s_name,'UL_CI_sim_events')
             obs_name=paste0(s_name,'obs_events')
             polygon(x=c(x_start-bin_width/2,x_start+bin_width/2,x_start+bin_width/2,x_start-bin_width/2),y=c(d[i,sim_LL_CI_name],d[i,sim_LL_CI_name],d[i,sim_UL_CI_name],d[i,sim_UL_CI_name]),col=col_scheme[k+2],border=fillcol)
             lines(x=c(x_start-bin_width/2,x_start+bin_width/2),y=c(d[i,obs_name],d[i,obs_name]),lwd=2,col=obs_lines) 
             lines(x=c(x_start-bin_width/2,x_start+bin_width/2),y=c(d[i,sim_med_name],d[i,sim_med_name]),lwd=1,col=sim_med) 
             # add vertical line indicating stratify separation
             if(add_vert==TRUE){
                xv=x_start+bin_width/2+strat_spacing/2
                next_sim_LL_CI_name=paste0(next_s_name,'LL_CI_sim_events')
                next_sim_UL_CI_name=paste0(next_s_name,'UL_CI_sim_events')
                yv1=min(d[i,sim_LL_CI_name],d[i,next_sim_LL_CI_name])
                yv2=max(d[i,sim_UL_CI_name],d[i,next_sim_UL_CI_name])
                lines(x=c(xv,xv),y=c(yv1,yv2),lty=3)
             }
           }
       }
     } 
     
     if(legend_position!=FALSE){
        if(stratify==FALSE){
            legend(legend_position,c('Observed event count','Median of simulated event counts',
                                     paste0(VPCdata$confidence_interval,'% Confidence interval of simulated events')),
                                     pch=c(NA,NA,15),pt.cex=c(NA,NA,2),lty=c(1,1,NA),lwd=c(2,1,NA),cex=legend_size,col=c(col_scheme[1],col_scheme[2],col_scheme[3]))
        }else{
          txt=c('Observed event count','Median of simulated event counts')
          for(m in 1:VPCdata$num_strat_levels){
            txt=c(txt,paste0(VPCdata$confidence_interval,'% Confidence interval of simulated events: ',VPCdata$stratification,'=',VPCdata$strat_values[m]))
          }
          col_vals=col_scheme[3:(2+VPCdata$num_strat_levels)]
          legend(legend_position,txt,
                 pch=c(NA,NA,rep(15,VPCdata$num_strat_levels)),pt.cex=c(NA,NA,rep(2,VPCdata$num_strat_levels)),lty=c(1,1,rep(NA,VPCdata$num_strat_levels)),lwd=c(2,1,rep(NA,VPCdata$num_strat_levels)),cex=legend_size,col=c(col_scheme[1],col_scheme[2],col_vals))
        }
     }
  }  # end of if(type%in%c('cumulative'))
  
  if(type%in%c('discrete')){
     d=VPCdata$values_noncumulative
     
     n_bins=nrow(d)
     if(stratify==TRUE){
        d$obs_events=NULL
        d$median_sim_events=NULL
        d$LL_CI_sim_events=NULL
        d$UL_CI_sim_events=NULL
     }
     plot(x=xrange,y=c(0,max(d[,-1])),type='n',xlab=x_label,ylab='Event count',...)
     
     if(stratify==FALSE){
        for(i in 1:n_bins){
          if (i==1){
            polygon(x=c(d$time[i]/2-bin_width/2,d$time[i]/2+bin_width/2,d$time[i]/2+bin_width/2,d$time[i]/2-bin_width/2),y=c(d$LL_CI_sim_events[i],d$LL_CI_sim_events[i],d$UL_CI_sim_events[i],d$UL_CI_sim_events[i]),col=fillcol,border=fillcol)
            lines(x=c(d$time[i]/2-bin_width/2,d$time[i]/2+bin_width/2),y=c(d$obs_events[i],d$obs_events[i]),lwd=2,col=obs_lines) 
            lines(x=c(d$time[i]/2-bin_width/2,d$time[i]/2+bin_width/2),y=c(d$median_sim_events[i],d$median_sim_events[i]),lwd=1,col=sim_med)
          } else {
            polygon(x=c((d$time[i]+d$time[i-1])/2-bin_width/2,(d$time[i]+d$time[i-1])/2+bin_width/2,(d$time[i]+d$time[i-1])/2+bin_width/2,(d$time[i]+d$time[i-1])/2-bin_width/2),y=c(d$LL_CI_sim_events[i],d$LL_CI_sim_events[i],d$UL_CI_sim_events[i],d$UL_CI_sim_events[i]),col=fillcol,border=fillcol)
            lines(x=c((d$time[i]+d$time[i-1])/2-bin_width/2,(d$time[i]+d$time[i-1])/2+bin_width/2),y=c(d$obs_events[i],d$obs_events[i]),lwd=2,col=obs_lines) 
            lines(x=c((d$time[i]+d$time[i-1])/2-bin_width/2,(d$time[i]+d$time[i-1])/2+bin_width/2),y=c(d$median_sim_events[i],d$median_sim_events[i]),lwd=1,col=sim_med) 
          }
          
        }
       
        abline(v=c(0,d$time),lty=2,col="darkgreen")
       
     }
     
     if(stratify==TRUE){
        for(k in 1:VPCdata$num_strat_levels){
           s_level=VPCdata$strat_values[k]
           s_name=paste0(s_level,'_')
           
           if(k>1){
              start_offset=start_offset+bin_width+strat_spacing
           }
           add_vert=FALSE
           if(k<VPCdata$num_strat_levels & add_strat_lines==TRUE){
              add_vert=TRUE
              next_s_level=VPCdata$strat_values[k+1]
              next_s_name=paste0(next_s_level,'_')
           }
           for(i in 1:n_bins){
             if (i==1){
               ctime=d$time[i]/2
             } else {
               ctime=(d$time[i]+d$time[i-1])/2
             }
              
              x_start<<-ctime+start_offset
              sim_med_name=paste0(s_name,'median_sim_events')
              sim_LL_CI_name=paste0(s_name,'LL_CI_sim_events')
              sim_UL_CI_name=paste0(s_name,'UL_CI_sim_events')
              obs_name=paste0(s_name,'obs_events')
              polygon(x=c(x_start-bin_width/2,x_start+bin_width/2,x_start+bin_width/2,x_start-bin_width/2),y=c(d[i,sim_LL_CI_name],d[i,sim_LL_CI_name],d[i,sim_UL_CI_name],d[i,sim_UL_CI_name]),col=col_scheme[k+2],border=fillcol)
              lines(x=c(x_start-bin_width/2,x_start+bin_width/2),y=c(d[i,obs_name],d[i,obs_name]),lwd=2,col=obs_lines) 
              lines(x=c(x_start-bin_width/2,x_start+bin_width/2),y=c(d[i,sim_med_name],d[i,sim_med_name]),lwd=1,col=sim_med) 
              # add vertical line indicating stratify separation
              if(add_vert==TRUE){
                 xv=x_start+bin_width/2+strat_spacing/2
                 next_sim_LL_CI_name=paste0(next_s_name,'LL_CI_sim_events')
                 next_sim_UL_CI_name=paste0(next_s_name,'UL_CI_sim_events')
                 yv1=min(d[i,sim_LL_CI_name],d[i,next_sim_LL_CI_name])
                 yv2=max(d[i,sim_UL_CI_name],d[i,next_sim_UL_CI_name])
                 lines(x=c(xv,xv),y=c(yv1,yv2),lty=3)
              }
           }
        }
       
       abline(v=c(0,d$time),lty=2,col="darkgreen")
       
     } 
     
     if(type=='discrete'){
        if(legend_position!=FALSE){
           if(stratify==FALSE){
              legend(legend_position,c('Observed event count','Median of simulated event counts',
                                       paste0(VPCdata$confidence_interval,'% Confidence interval of simulated events')),
                     pch=c(NA,NA,15),pt.cex=c(NA,NA,2),lty=c(1,1,NA),lwd=c(2,1,NA),cex=legend_size,col=c(col_scheme[1],col_scheme[2],col_scheme[3]))
           }else{
              txt=c('Observed event count','Median of simulated event counts')
              for(m in 1:VPCdata$num_strat_levels){
                 txt=c(txt,paste0(VPCdata$confidence_interval,'% Confidence interval of simulated events: ',VPCdata$stratification,'=',VPCdata$strat_values[m]))
              }
              col_vals=col_scheme[3:(2+VPCdata$num_strat_levels)]
              legend(legend_position,txt,
                     pch=c(NA,NA,rep(15,VPCdata$num_strat_levels)),pt.cex=c(NA,NA,rep(2,VPCdata$num_strat_levels)),lty=c(1,1,rep(NA,VPCdata$num_strat_levels)),lwd=c(2,1,rep(NA,VPCdata$num_strat_levels)),cex=legend_size,col=c(col_scheme[1],col_scheme[2],col_vals))
           }
        }
     }
  }  # end of if(type%in%c('discrete','joint'))
}
##############################################################
####### END OF FUNCTION plotEventVPC #########################
##############################################################


##############################################################
### Function to adjust simulated event times to the next scheduled observation appointment for each patient
###
### Function Inputs:
### event_data_frame = a data frame of event times, event types and (optionally) patient IDs
### time_col = name of column containing the simulated event time
### event_col = name of column containing event types (1=event, 0=dropout)
### patient_num_col = name of column containing patient number
### appointment_times = a list of vectors of appointment times for each patient
### time_for_events_beyond_last_appointment = if a simulated event is at a time beyond the final appointment then assign the adjusted event time to be time_for_events_beyond_last_appointment
### adjust_censorings (T/F). Determines if censoring events should also have event time shifted. For example, this would be set to T if censoring from a tumour progression model was also interval censored because events were censored at times of tumour assessments
###
### Function Outputs:
### the function returns the event_data_frame that was input, but with additional column added called adjusted_event_time
##############################################################
adjust_events_to_appointment_times=function(event_data_frame,time_col='event_time',event_col='event',
            patient_ID_col='patient_ID',appointment_times,time_for_events_beyond_last_appointment,adjust_censorings=FALSE){

  event_data_frame$adjusted_event_time=-9
  for(i in 1:nrow(event_data_frame)){
    new_event_time=-9
    current_time=event_data_frame[i,time_col]
    current_event_type=event_data_frame[i,event_col]
    current_patient=as.character(event_data_frame[i,patient_ID_col])
    if(current_event_type==1 | adjust_censorings==TRUE){
      # only adjust events. don't adjust censorings
      current_appointment_times=appointment_times[[current_patient]]
      new_event_time=current_appointment_times[current_appointment_times>current_time][1]
      if(is.na(new_event_time)==TRUE){
        new_event_time=time_for_events_beyond_last_appointment
      }
    } else {
      new_event_time=current_time
    }
    event_data_frame[i,'adjusted_event_time']=new_event_time
  }
  return(event_data_frame)
}
##############################################################
####### END OF FUNCTION adjust_events_to_appointment_times ###
##############################################################



##############################################################
### Function to adjust repeated simulated event times to the next scheduled observation appointment for each patient
###
### Function Inputs:
### event_list = a list of data frames of event times, event types and patient IDs
### num_sim_repeats = number of simulation repeats
### time_col = name of column containing the simulated event time
### event_col = name of column containing event types (1=event, 0=dropout)
### patient_ID_col = name of column containing patient number
### appointment_times = either a list (of length num_patients) of vectors of appointment times for each patient, when we want to apply same appointment times for each simulation repeat
###                     or it can be a list (of length num_sim_repeats) which contains num_patients sub-lists, each of which contains vectors of appointment times for each patient, when we want to apply different appointment times for each simulation repeat
### time_for_events_beyond_last_appointment = if a simulated event is at a time beyond the final appointment then assign the adjusted event time to be time_for_events_beyond_last_appointment
### adjust_censorings (T/F). Determins if censoring events should also be shifted. For example, this would be set to T if censoring from a tumour progression model was also interval censoed because events were censored at times of tumour assessments
###
### Function Outputs:
### the function returns the event_list that was input, but with each constituent data frame containing an additional column added called adjusted_event_time
##############################################################
adjust_repeat_events_to_appointment_times=function(event_list,num_sim_repeats,time_col='event_time',event_col='event',
             patient_ID_col='patient_ID',appointment_times,time_for_events_beyond_last_appointment,adjust_censorings=FALSE){
  
  adjusted_event_list=event_list
  print('Adjusting events to times of assessments...')
  for(i in 1:num_sim_repeats){
    #print(i)
    # here we figure out if the appointment_times is a list of vectors (same times for each simulation repeat) or a list of lists of vectors (different times for each simulation repeat)
    if(length(appointment_times[[1]][[1]])<2){
      appointment_times_n=appointment_times  
    } else {
      appointment_times_n=appointment_times[[i]]
    }

    # call the adjust_events_to_appointment_times function
    event_data_frame=event_list$events[[i]]
    result=adjust_events_to_appointment_times(event_data_frame,time_col,event_col,patient_ID_col,appointment_times_n,time_for_events_beyond_last_appointment,adjust_censorings)
    adjusted_event_list$events[[i]]=result
  }
  return(adjusted_event_list)
}
#####################################################################
####### END OF FUNCTION adjust_repeat_events_to_appointment_times ###
#####################################################################

##############################################################
### Function to calculate confidence interval as central X% range of many simulated Kaplan-Meier curves
###
### Function Inputs:
### event_list = a list of data frames (one for each simulation repeat) of simulated event times and event types
###              it will have been generated using the repeated_sim_events function
### time_col = name of column (within the data frames in event_list) containing event times
### event_col = name of column (within the data frames in event_list) containing event type (1=event,0=censored)
### time_points = a vector of times at which each of the many KM curves will be interpolated. Must begin at time=zero
### conf_range_as_percent = the central percentage of interpolated KM curves that form the confidence bands
###
### Function Outputs:
### A list of 2 data frames:
### the first is a data from of interpolated survival curves for each simulation repeat
### the second is a data frame containing time, lower confidence band, upper confidence band, median
##############################################################
calc_KM_conf_bands=function(event_list,time_col='event_time',event_col='event',time_points,conf_range_as_percent=95){
   
   num_sim_repeats=length(event_list$events)
   
   print('Calculating confidence bands around KM curve...')
   
   # calculate a KM curve for each repeated simulation  
   KM_list=vector(mode='list',length=num_sim_repeats)
   for(i in 1:num_sim_repeats){
      E_times=event_list$events[[i]][,time_col]
      E_types=event_list$events[[i]][,event_col]
      s1=Surv(E_times,E_types)
      km1=survfit(s1~1,conf.int=0.95)
      KM_list[[i]]=km1
   }
   
   # create list of event times from each simulation repeat
   sim_times=vector(mode='list',length=num_sim_repeats)
   for (k in 1:num_sim_repeats) {
      sim_times[[k]]=KM_list[[k]]$time    
   }
   
   # create list of survival probabilities from each simulation repeat
   sim_probs=vector(mode='list',length=num_sim_repeats)
   for (k in 1:num_sim_repeats) {
      sim_probs[[k]]=KM_list[[k]]$surv   
   }  
   
   # create a vector of times for our interpolations from the simulated Kaplan-Meier curves
   t1=time_points
   
   # create a data frame to contain these times plus interpolated probabilities from each simulation
   # first column is the time of interpolation of events
   # each additional column contains the interpolated Kaplan-Meier probability from each simulation (initially set to -9)
   x=rep(-9,length(t1))
   sim_results=data.frame(t1)
   for (q in 1:num_sim_repeats){
      sim_results=data.frame(sim_results,x) 
   }
   
   # now loop through each simulation and at each time in vector t1, interpolate 
   # and put results into data frame sim_results
   # the survival probability
   for (i in 1:num_sim_repeats) {
      prob_vals=NULL
      prob_vals[1]=1      # when time=0 probability=1
      times_to_search=sim_times[[i]]
      for (k in 2:length(t1)){
         current_time=t1[k]
         # at this current_time we need to search the relevant time column in sim_times
         # for the first time value that is less then or equal to current time, then look
         # up the probability at that time in sim_probs
         times1=times_to_search[times_to_search<=current_time]
         index=length(times1)
         if (index==0){
            current_prob=1
         } else {
            current_prob=sim_probs[[i]][index]
         }
         prob_vals[k]=current_prob  
      }
      sim_results[,i+1]=prob_vals
   }
   
   # for each time point in vector t1 we now calculate the
   # percentiles and save these values
   # in a new data frame called final_results
   
   # create the final_results data frame and initial fill the confidence intervals with -9
   final_results=data.frame(t1)
   final_results$lower_conf=rep(-9,length(t1))
   final_results$upper_conf=rep(-9,length(t1))
   final_results$median=rep(-9,length(t1))
   
   # now calculate the central x% range at each time point of vector t1
   upper_val=conf_range_as_percent/100 + (1 - conf_range_as_percent/100)/2
   lower_val=(1 - conf_range_as_percent/100)/2
   for (i in 1:length(t1)){
      q=sim_results[i,2:(num_sim_repeats+1)]
      c_upper=quantile(q,probs=upper_val)
      c_lower=quantile(q,probs=lower_val)
      c_med=median(as.numeric(q)) #   quantile(q,probs=0.5)
      final_results[i,'upper_conf']=c_upper
      final_results[i,'lower_conf']=c_lower
      final_results[i,'median']=c_med 
   }
   
   res=list(interpolated_survival_curves=sim_results,confidence_bands=final_results)
   print('Confidence bands have been calculated')
   return(res)
   
}
##############################################################
####### END OF FUNCTION calc_KM_conf_bands ###################
############################################################## 


###############################################################################
### FUNCTION FOR ADDING SIMULATION CONFIDENCE BANDS TO AN EXISTING KM PLOT ####
###############################################################################
### LIST OF FUNCTION ARGUMENTS
#
#  conf_band_data_frame     # data frame containing the calculated confidence bands - this will be the output of the calc_KM_conf_bands function
#  band_colour              # colour for the confidence bands, set to NA to exclude
#  median_colour            # colour for the median curve, set to NA to exclude
#  band_line_width          # line width for the confidence bands (typically 1 or 2)
#  median_line_width        # line width for the median curve (typically 1 or 2)
#  band_density             # density for the confidence region (typically 10-20)
#  band_angle               # angle for the shading of the confidence bands (typicaly -45 to +45)
###############################################################################
add_KM_conf_bands=function(conf_band_data_frame,band_colour='gray',median_colour='red',
                           band_line_width=1,median_line_width=2,band_density=100,band_angle=0){
   
   if(is.na(band_colour)==FALSE){
      w1=c(conf_band_data_frame[,1],rev(conf_band_data_frame[,1]))
      w2=c(conf_band_data_frame$upper_conf,rev(conf_band_data_frame$lower_conf))
      polygon(x=w1,y=w2,col=band_colour,border=band_colour,density=band_density,angle=band_angle,lwd=band_line_width)
   }
   if(is.na(median_colour)==FALSE){  
      points(x=conf_band_data_frame[,1],y=conf_band_data_frame$median, col=median_colour,type='l',lwd=median_line_width)
   }
}
###############################################################################
###### END OF FUNCTION add_KM_conf_bands ######################################
###############################################################################

