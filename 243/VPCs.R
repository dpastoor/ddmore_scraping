###################################################################################################
###################################################################################################
###################################################################################################
###
### BAST TTE modelling 
### VPC generation example script
###
###################################################################################################
###################################################################################################
###################################################################################################

rm(list=ls(all=TRUE))          # Reset memory

# Set working directory and load BAST survival function library 
setwd("//BAST_SERVER/BAST_Shared_Files/BAST General/BAST Publications/DDMoRe Repository/TTE of competing risks/Final documents/R scripts/VPCs/")    # Set your own work directory
source('..\\BAST_surv_functions.R')
set.seed(28674)

# set the ID column name in the data files we will work with
ID_col='ID'

# load the event data
all_data=read.csv("..\\..\\Data\\event_data.csv")

# load the times of each patient's assessment appointments. These are required for simulations of interval-censored events, where simulated
# event times need to be shifted to the time of the next assessment apppointment.
v_data=read.table("..\\..\\Data\\assessment_times.csv",header=TRUE,sep=",",stringsAsFactors=FALSE) 

#######################################################################################################################
### DEFINE MODELS ###
pars_EV1=c(HazardFunction='exponential',lambda=2.8/1000,alpha=NA,CovMod=1) 
pars_EV2=c(HazardFunction='gompertz',lambda=3.28/1000,alpha=4.05/1000,CovMod=1)    
pars_COMPEV1=c(HazardFunction='log_normal',lambda=1.42,alpha=984,CovMod=1) 
pars_COMPEV2=c(HazardFunction='log_normal',lambda=1.08,alpha=167,CovMod=0) 
#######################################################################################################################

#######################################################################################################################
#######################################################################################################################
### LOAD DATA

# EV1 data
EV1_data = all_data[all_data$DVID %in% c(0,1),]
EV1_data = EV1_data[order(EV1_data[,ID_col],EV1_data[,'TIME']),]  # make sure data is ordered according to the patient identifier we are using
EV1_data = EV1_data[EV1_data$DV %in% c(-1,1),]
EV1_data[EV1_data$DV==-1,'DV']=0

# EV2 data
EV2_data = all_data[all_data$DVID %in% c(0,2),]
EV2_data = EV2_data[order(EV2_data[,ID_col],EV2_data[,'TIME']),]  # make sure data is ordered according to the patient identifier we are using
EV2_data = EV2_data[EV2_data$DV %in% c(-1,3),]
EV2_data[EV2_data$DV==-1,'DV']=0
EV2_data[EV2_data$DV==3,'DV']=1

# COMPEV1 data
COMPEV1_data = all_data[all_data$DVID %in% c(0,3),]
COMPEV1_data = COMPEV1_data[order(COMPEV1_data[,ID_col],COMPEV1_data[,'TIME']),]  # make sure data is ordered according to the patient identifier we are using
COMPEV1_data = COMPEV1_data[COMPEV1_data$DV %in% c(-1,3),]
COMPEV1_data[COMPEV1_data$DV==-1,'DV']=0
COMPEV1_data[COMPEV1_data$DV==3,'DV']=1

# COMPEV2 data
COMPEV2_data = all_data[all_data$DVID %in% c(0,4),]
COMPEV2_data = COMPEV2_data[order(COMPEV2_data[,ID_col],COMPEV2_data[,'TIME']),]  # make sure data is ordered according to the patient identifier we are using
COMPEV2_data = COMPEV2_data[COMPEV2_data$DV %in% c(-1,1),]
COMPEV2_data[COMPEV2_data$DV==-1,'DV']=0
#######################################################################################################################

#######################################################################################################################
### PRODUCE KM PLOTS FOR THE FOUR EVENT TYPES 
time_range=c(0,1000)      # set time range for x-axis of KM plots
par(mfrow=c(2,2))
plot_KM(EV1_data,time_col='TIME',x_range=time_range,main='Base model for EV1_data')
surv_EV1=IntHaz(parms=c(HazardFunction='exponential',lambda=3.21/1000,alpha=NA,CovMod=0),maxtime=time_range[2],timeint=1)   # calculate survival curve for base model using numerical integration of hazard function
points(x=surv_EV1$time,y=surv_EV1$SURV,type='l',col='red',lwd=2)                                                             # add the calculated base model survival curve to plot

plot_KM(EV2_data,time_col='TIME',x_range=time_range,main='Base model for EV2_data')
surv_EV2=IntHaz(parms=c(HazardFunction='gompertz',lambda=3.54/1000,alpha=3.32/1000,CovMod=0),maxtime=time_range[2],timeint=1)   
points(x=surv_EV2$time,y=surv_EV2$SURV,type='l',col='red',lwd=2)  

plot_KM(COMPEV1_data,time_col='TIME',x_range=time_range,main='Base model for COMPEV1_data')
surv_CEV1=IntHaz(parms=c(HazardFunction='log_normal',lambda=1.29,alpha=620,CovMod=0),maxtime=time_range[2],timeint=1)   
points(x=surv_CEV1$time,y=surv_CEV1$SURV,type='l',col='red',lwd=2) 

plot_KM(COMPEV2_data,time_col='TIME',x_range=time_range,main='Base model for COMPEV2_data')
surv_CEV2=IntHaz(parms=pars_COMPEV2,maxtime=time_range[2],timeint=1)   
points(x=surv_CEV2$time,y=surv_CEV2$SURV,type='l',col='red',lwd=2) 
#######################################################################################################################

#######################################################################################################################
### SETUP PATIENT INFO DATAFRAME
# we will need a data frame of IDs and associated database lock times
# Also include covariate categories to be used in the stratified VPCs later
inf_patient=all_data[!duplicated(all_data$ID),c('ID','DBLOCK')]
names(inf_patient)=c('ID','DB_LOCK_TIME')
IDs=unique(all_data[,ID_col])
n_IDs=length(IDs)
for (i in 1:n_IDs){
  inf_patient$AGE_CAT[inf_patient$ID==IDs[i]]=all_data$AGE_CAT[all_data$EVID==2 & all_data$ID==IDs[i]]
  inf_patient$NEUT_CAT[inf_patient$ID==IDs[i]]=all_data$NEUT_CAT[all_data$EVID==2 & all_data$ID==IDs[i]]
  inf_patient$AUC_CAT[inf_patient$ID==IDs[i]]=all_data$AUC_CAT[all_data$EVID==2 & all_data$ID==IDs[i]]
}


#######################################################################################################################
#######################################################################################################################
### CALCULATION OF COVARIATE INFLUENCE ################################################################################
# In this example the covariate influences do not vary with time. However, the function library can deal with time-varying
# covariates. Hence, if the values of NEUT, AGE and AUC within all_data did vary with time, the code below will still work
# and will incorporate the time varying covariate influences during numerical intergration of the hazard functions.


ls_cov_EV1=vector(mode='list',length=n_IDs)     # we will define the influence of NEUT and AGE on hazard of EV1 for each patient and save to this object
ls_cov_EV2=vector(mode='list',length=n_IDs)     # we will define the influence of AUC on hazard of EV2 for each patient and save to this object
ls_cov_COMPEV1=vector(mode='list',length=n_IDs) # we will define the influence of AGE on hazard of COMPEV1 for each patient and save to this object
ls_t=vector(mode='list',length=n_IDs)           # the corresponding times will be saved to this object
names(ls_t)=as.character(IDs)
names(ls_cov_EV1)=as.character(IDs)
names(ls_cov_EV2)=as.character(IDs)
names(ls_cov_COMPEV1)=as.character(IDs)

# loop through each patient and calculate covariate influence(s) and save to the list objects: ls_cov_EV1, ls_cov_EV2 and ls_cov_COMPEV1
for (i in 1:n_IDs) {
   id_dat=all_data[all_data[,ID_col]==IDs[i],]                                         # get data for current patient
   ls_cov_EV1[[IDs[i]]]=exp(-0.000156*(id_dat$NEUT-4133))*exp(0.032*(id_dat$AGE-55))  # influence of NEUT and AGE on hazard of EV1
   ls_cov_EV2[[IDs[i]]]=exp(0.000309*(id_dat$AUC-3065.5))                              # influence of AUC on hazard of EV2
   ls_cov_COMPEV1[[IDs[i]]]=exp(0.0509*(id_dat$AGE-55))                                 # influence of AGE on hazard of COMPEV1 
   ls_t[[IDs[i]]]=id_dat$TIME                                                          # time values
}
#######################################################################################################################

#############################################################################################################################################
### RUN SIMULATIONS
# calculate individual hazard and survival curves for each patient
surv_EV1=sim_events(Patient_Info=inf_patient,Int_Time_Col='DB_LOCK_TIME',Parameters=pars_EV1,ls_times=ls_t,ls_cov=ls_cov_EV1)
surv_EV2=sim_events(Patient_Info=inf_patient,Int_Time_Col='DB_LOCK_TIME',Parameters=pars_EV2,ls_times=ls_t,ls_cov=ls_cov_EV2)
surv_COMPEV1=sim_events(Patient_Info=inf_patient,Int_Time_Col='DB_LOCK_TIME',Parameters=pars_COMPEV1,ls_times=ls_t,ls_cov=ls_cov_COMPEV1)
surv_COMPEV2=sim_events(Patient_Info=inf_patient,Int_Time_Col='DB_LOCK_TIME',Parameters=pars_COMPEV2)

# simulate 1000 sets of event times for each patient.
event_EV1=repeated_sim_events(Patient_Info=inf_patient,num_sim_repeats=1000,survival_curves=surv_EV1$surv_haz)
event_EV2=repeated_sim_events(Patient_Info=inf_patient,num_sim_repeats=1000,survival_curves=surv_EV2$surv_haz)
event_COMPEV1=repeated_sim_events(Patient_Info=inf_patient,num_sim_repeats=1000,survival_curves=surv_COMPEV1$surv_haz)
event_COMPEV2=repeated_sim_events(Patient_Info=inf_patient,num_sim_repeats=1000,survival_curves=surv_COMPEV2$surv_haz)
#############################################################################################################################################

#############################################################################################################################################
### BASIC QC. ENSURE THAT KM PLOTS OF FIRST 40 EVENT SIMULATIONS LOOKS SIMILAR TO OBSERVED EVENTS
par(mfrow=c(2,2))
plot_KM(EV1_data,time_col='TIME',main='EV1 model',x_range=time_range)
for (i in 1:40){
  plot_KM(event_EV1$events[[i]],x_range=time_range,time_col='event_time',event_col='event',KM_colour='red',mark_censored_events=FALSE,add_to_existing_plot=TRUE,conf_interval_pct=FALSE,line_width=1)
}
plot_KM(EV1_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
plot_KM(EV2_data,time_col='TIME',main='EV2 model',x_range=time_range)
for (i in 1:40){
  plot_KM(event_EV2$events[[i]],x_range=time_range,time_col='event_time',event_col='event',KM_colour='red',mark_censored_events=FALSE,add_to_existing_plot=TRUE,conf_interval_pct=FALSE,line_width=1)
}
plot_KM(EV2_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
plot_KM(COMPEV1_data,time_col='TIME',main='COMPEV1 model',x_range=time_range)
for (i in 1:40){
  plot_KM(event_COMPEV1$events[[i]],x_range=time_range,time_col='event_time',event_col='event',KM_colour='red',mark_censored_events=FALSE,add_to_existing_plot=TRUE,conf_interval_pct=FALSE,line_width=1)
}
plot_KM(COMPEV1_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
plot_KM(COMPEV2_data,time_col='TIME',main='COMPEV2 model',x_range=time_range)
for (i in 1:40){
  plot_KM(event_COMPEV2$events[[i]],x_range=time_range,time_col='event_time',event_col='event',KM_colour='red',mark_censored_events=FALSE,add_to_existing_plot=TRUE,conf_interval_pct=FALSE,line_width=1)
}
plot_KM(COMPEV2_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
#############################################################################################################################################

#############################################################################################################################################
#### SHIFT SIMULATED EVENTS TO NEXT TIME OF ASSESSMENT
# EV2 and COMPEV1 are of interval-censored event type, and the simulated event times need to be shifted to the next assessment visit time.
# Build the assessment visit times for each patient into an appt_time_list object
appt_time_list=vector(mode='list',length=length(IDs))
for(i in 1:length(IDs)){
  # extract the attended visits for patient i  
  current_data=v_data[v_data$ID==IDs[i],'TIME']
  # Need to add planned assessment times up to time of database lock for this patient, assume every 60 days
  DB_lock_time=inf_patient$DB_LOCK_TIME[inf_patient$ID==IDs[i]]
  if(length(current_data)>0){
    extra_times=seq(tail(current_data,1),DB_lock_time,60)
  } else {
    extra_times=seq(0,DB_lock_time,60)
  }
  appt_time_list[[i]]=unique(c(current_data,extra_times))
}
names(appt_time_list)=as.character(IDs)

# adjust simulated event times to next scheduled assessment.
event_EV2=adjust_repeat_events_to_appointment_times(event_list=event_EV2,num_sim_repeats=1000,appointment_times=appt_time_list,time_for_events_beyond_last_appointment=10000)
event_COMPEV1=adjust_repeat_events_to_appointment_times(event_list=event_COMPEV1,num_sim_repeats=1000,appointment_times=appt_time_list,time_for_events_beyond_last_appointment=10000)

# now visualize the effect of shifting event times for EV2 and COMPEV1
par(mfrow=c(2,2))
plot_KM(EV2_data,time_col='TIME',main='EV2 model',x_range=time_range)
for (i in 1:40){
   plot_KM(event_EV2$events[[i]],x_range=time_range,time_col='event_time',event_col='event',KM_colour='red',mark_censored_events=FALSE,add_to_existing_plot=TRUE,conf_interval_pct=FALSE,line_width=1)
}
plot_KM(EV2_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
plot_KM(EV2_data,time_col='TIME',main='EV2 model (adjusted sim events)',x_range=time_range)
for (i in 1:40){
   plot_KM(event_EV2$events[[i]],x_range=time_range,time_col='adjusted_event_time',event_col='event',KM_colour='red',mark_censored_events=FALSE,add_to_existing_plot=TRUE,conf_interval_pct=FALSE,line_width=1)
}
plot_KM(EV2_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
plot_KM(COMPEV1_data,time_col='TIME',main='COMPEV1 model',x_range=time_range)
for (i in 1:40){
   plot_KM(event_COMPEV1$events[[i]],x_range=time_range,time_col='event_time',event_col='event',KM_colour='red',mark_censored_events=FALSE,add_to_existing_plot=TRUE,conf_interval_pct=FALSE,line_width=1)
}
plot_KM(COMPEV1_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
plot_KM(COMPEV1_data,time_col='TIME',main='COMPEV1 model (adjusted sim events)',x_range=time_range)
for (i in 1:40){
   plot_KM(event_COMPEV1$events[[i]],x_range=time_range,time_col='adjusted_event_time',event_col='event',KM_colour='red',mark_censored_events=FALSE,add_to_existing_plot=TRUE,conf_interval_pct=FALSE,line_width=1)
}
plot_KM(COMPEV1_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
#############################################################################################################################################

#############################################################################################################################################
### RUN COMPETING EVENT ANALYSIS
# define a list object that contains, for each event type, all simulated events that compete with it. Names used in the list definition will be used in the function output
comp1=list(COMPEV1=event_COMPEV1,COMPEV2=event_COMPEV2)   # For use with both EV1 and EV2: COMPEV1 and COMPEV2 compete
comp2=list(COMPEV1=event_COMPEV1)                         # For use with COMPEV2: COMPEV1 competes
comp3=list(COMPEV2=event_COMPEV2)                         # For use with COMPEV1: COMPEV2 competes

# Run the competing events analysis. If a competing event is simulated before the event of interest then the event of interest becomes censored at the time of the competing event.
final_EV1=CompeteEvents(event_name='EV1',events=event_EV1,competing_events=comp1)
final_EV2=CompeteEvents(event_name='EV2',events=event_EV2,competing_events=comp1)
final_COMPEV1=CompeteEvents(event_name='COMPEV1',events=event_COMPEV1,competing_events=comp3)
final_COMPEV2=CompeteEvents(event_name='COMPEV2',events=event_COMPEV2,competing_events=comp2)
#############################################################################################################################################

#############################################################################################################################################
### CALCULATE  EVENT VPCs (this can be done before or after running CompeteEvents(). If done after running CompeteEvents() then toggle compete=TRUE/FALSE to see effect of competing events)
VPC_EV2=calcEventVPC(Patient_Info=inf_patient,obs_events=EV2_data,obs_event_cols=c('TIME','DV'),sim_events=final_EV2,compete=TRUE,bins=c(10,60,100,200,300,400),conf_int=90,stratify="AUC_CAT")
VPC_COMPEV1=calcEventVPC(Patient_Info=inf_patient,obs_events=COMPEV1_data,obs_event_cols=c('TIME','DV'),sim_events=final_COMPEV1,compete=TRUE,bins=c(10,60,100,200,300,400),conf_int=90,stratify="AGE_CAT")
VPC_COMPEV2=calcEventVPC(Patient_Info=inf_patient,obs_events=COMPEV2_data,obs_event_cols=c('TIME','DV'),sim_events=final_COMPEV2,compete=TRUE,bins=c(10,60,100,200,300,400),conf_int=90,stratify=FALSE)

# First split VPC of all patients by NEUT_CAT
# Then, split again by AGE_CAT
VPC_EV1_1=calcEventVPC(Patient_Info=inf_patient,obs_events=EV1_data,obs_event_cols=c('TIME','DV'),sim_events=final_EV1,compete=TRUE,bins=c(10,60,100,200,300,400),conf_int=90,stratify="NEUT_CAT")
VPC_EV1_2=calcEventVPC(Patient_Info=inf_patient,obs_events=EV1_data,obs_event_cols=c('TIME','DV'),sim_events=final_EV1,compete=TRUE,bins=c(10,60,100,200,300,400),conf_int=90,stratify="AGE_CAT")
#############################################################################################################################################

#############################################################################################################################################
# SAVE VPC PLOTS - cumulative
png(file='VPC_EV1_1.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_EV1_1,display_percent=TRUE,type='cumulative',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
png(file='VPC_EV1_2.png',width=800,height=650) 
plotEventVPC(VPCdata=VPC_EV1_2,display_percent=TRUE,type='cumulative',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()

png(file='VPC_EV2.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_EV2,display_percent=TRUE,type='cumulative',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
png(file='VPC_EV2_nostrat.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_EV2,display_percent=TRUE,type='cumulative',stratify=FALSE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()


png(file='VPC_COMPEV1.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_COMPEV1,display_percent=TRUE,type='cumulative',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
png(file='VPC_COMPEV2.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_COMPEV2,display_percent=TRUE,type='cumulative',stratify=FALSE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()

#############################################################################################################################################
# SAVE VPC PLOTS - discrete
png(file='VPC_EV1_1_dis.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_EV1_1,display_percent=TRUE,type='discrete',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
png(file='VPC_EV1_2_dis.png',width=800,height=650) 
plotEventVPC(VPCdata=VPC_EV1_2,display_percent=TRUE,type='discrete',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()

png(file='VPC_EV2_dis.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_EV2,display_percent=TRUE,type='discrete',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
png(file='VPC_EV2_nostrat_dis.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_EV2,display_percent=TRUE,type='discrete',stratify=FALSE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()


png(file='VPC_COMPEV1_dis.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_COMPEV1,display_percent=TRUE,type='discrete',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
png(file='VPC_COMPEV2_dis.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_COMPEV2,display_percent=TRUE,type='discrete',stratify=FALSE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
#############################################################################################################################################
  
#############################################################################################################################################
### NOW REPEAT VPCs BUT WITHOUT ANALYSIS OF COMPETING EVENTS
### CALCULATE  EVENT VPCs
VPC_EV2a=calcEventVPC(Patient_Info=inf_patient,obs_events=EV2_data,obs_event_cols=c('TIME','DV'),sim_events=final_EV2,compete=FALSE,bins=c(10,60,100,200,300,400),conf_int=90,stratify="AUC_CAT")
VPC_COMPEV1a=calcEventVPC(Patient_Info=inf_patient,obs_events=COMPEV1_data,obs_event_cols=c('TIME','DV'),sim_events=final_COMPEV1,compete=FALSE,bins=c(10,60,100,200,300,400),conf_int=90,stratify="AGE_CAT")
VPC_COMPEV2a=calcEventVPC(Patient_Info=inf_patient,obs_events=COMPEV2_data,obs_event_cols=c('TIME','DV'),sim_events=final_COMPEV2,compete=FALSE,bins=c(10,60,100,200,300,400),conf_int=90,stratify=FALSE)

# EV1, split by NEUT then AGE
VPC_EV1_1a=calcEventVPC(Patient_Info=inf_patient,obs_events=EV1_data,obs_event_cols=c('TIME','DV'),sim_events=final_EV1,compete=FALSE,bins=c(10,60,100,200,300,400),conf_int=90,stratify="NEUT_CAT")
VPC_EV1_2a=calcEventVPC(Patient_Info=inf_patient,obs_events=EV1_data,obs_event_cols=c('TIME','DV'),sim_events=final_EV1,compete=FALSE,bins=c(10,60,100,200,300,400),conf_int=90,stratify="AGE_CAT")
#############################################################################################################################################

#############################################################################################################################################
# PRINT VPC PLOTS
png(file='VPC_EV1_1_no_competition.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_EV1_1a,display_percent=TRUE,type='cumulative',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()

png(file='VPC_EV1_2_no_competition.png',width=800,height=650) 
plotEventVPC(VPCdata=VPC_EV1_2a,display_percent=TRUE,type='cumulative',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()

png(file='VPC_EV2_no_competition.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_EV2a,display_percent=TRUE,type='cumulative',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
png(file='VPC_COMPEV1_no_competition.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_COMPEV1a,display_percent=TRUE,type='cumulative',stratify=TRUE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
png(file='VPC_COMPEV2_no_competition.png',width=800,height=650)
plotEventVPC(VPCdata=VPC_COMPEV2a,display_percent=TRUE,type='cumulative',stratify=FALSE,bin_width=10,strat_spacing=3,legend_position="topright",xlim=c(0,400),ylim=c(0,100))
dev.off()
#############################################################################################################################################

#############################################################################################################################################
### NOW COMPARE EVENT VPCs WITH THE KAPLAM-MEIER VIEW (example using EV1)

# calculate unstratified event VPCs with and without consideration of competing events
VPC_EV1_no_strat_no_compete=calcEventVPC(Patient_Info=inf_patient,obs_events=EV1_data,obs_event_cols=c('TIME','DV'),sim_events=final_EV1,compete=FALSE,bins=c(10,60,100,200,300,400),conf_int=90,stratify=FALSE)
VPC_EV1_no_strat_compete=calcEventVPC(Patient_Info=inf_patient,obs_events=EV1_data,obs_event_cols=c('TIME','DV'),sim_events=final_EV1,compete=TRUE,bins=c(10,60,100,200,300,400),conf_int=90,stratify=FALSE)

# plot unstratified event VPCs with and without consideration of competing events.
# Note how there are fewer simulated events when competing events are considered and that this better matches the observed event counts.
png(file='z_event_VPC_EV1_no_strat_with_compete.png',width=800,height=650) 
plotEventVPC(VPCdata=VPC_EV1_no_strat_compete,display_percent=TRUE,type='cumulative',main='EV1 simulations with consideration of competing events',stratify=FALSE,bin_width=10,strat_spacing=3,legend_position="bottomright",xlim=c(0,400),ylim=c(0,100))
dev.off()
png(file='z_event_VPC_EV1_no_strat_without_compete.png',width=800,height=650) 
plotEventVPC(VPCdata=VPC_EV1_no_strat_no_compete,display_percent=TRUE,type='cumulative',main='EV1 simulations without consideration of competing events',stratify=FALSE,bin_width=10,strat_spacing=3,legend_position="bottomright",xlim=c(0,400),ylim=c(0,100))
dev.off()

# calculate confidence intervals around the 1000 simulated KM curves, with and without consideration of competing events
conf_KM_EV1_no_strat_no_compete=calc_KM_conf_bands(event_list=final_EV1,time_col='event_time',event_col='event',time_points=0:1000,conf_range_as_percent=90)
conf_KM_EV1_no_strat_compete=calc_KM_conf_bands(event_list=final_EV1,time_col='new_event_time',event_col='new_event',time_points=0:1000,conf_range_as_percent=90)

# plot observed and simulated KM curves, with and without consideration of competing events
### IMPORTANT. Note how the the median of the simulated KM curves is almost identical both with and without consideration of competing events.
### This is because the KM estimator seeks to determine the underlying event rate. When competing events are considered, some of the simulated 
### events are converted to censorings but the underlying event rate (produced by integration of the hazard function) is the same.
### when producing event simulations without consideration of competing events, the good fit of observed and simulated data in KM VPC can lead one 
### to the wrong conclusion that the simulated event counts produce a correct picture of expected event counts. Only with the consideration of 
### competing events, along with visual demonstration using the event VPC, will realistic event counts be generated and confirmed.
png(file='z_KM_VPC_EV1_no_strat_with_without_competition.png',width=1200,height=650) 
par(mfrow=c(1,2))
# KM plot where competing events have been analysed
plot_KM(EV1_data,time_col='TIME',main='EV1 simulations with consideration of competing events',x_range=time_range)
add_KM_conf_bands(conf_KM_EV1_no_strat_compete$confidence_bands)
add_to_existing_plot=TRUE
plot_KM(EV1_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
# KM plot where competing events have been ignored
plot_KM(EV1_data,time_col='TIME',main='EV1 simulations without consideration of competing events',x_range=time_range)
add_KM_conf_bands(conf_KM_EV1_no_strat_no_compete$confidence_bands)
add_to_existing_plot=TRUE
plot_KM(EV1_data,time_col='TIME',x_range=time_range,add_to_existing_plot=TRUE)
dev.off()
#############################################################################################################################################

save.image(file='R_objects.RData')

