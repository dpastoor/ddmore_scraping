#' cipro models for course infectious diseases (2016)
#' =================================================

#' TASK2 start (Set files and paths, data exploration)
#' ===========
pName     <- "HO_Bacteria_PKPD"
mDir      <- "models"
datafile  <- "cipro202.csv"
cipromdl  <- "cipro202.mdl"
modelwd   <- file.path(.MDLIDE_WORKSPACE_PATH, pName, mDir)

#' Set root directory
setwd(modelwd)
#opts_knit$set(root.dir = modelwd)      #if you want to run knitr and create a summary
getwd()

# PLOT data task 2, or USE excel, your own R-script, xpose etc or provided script below

#' Summarize data
data1  <- read.csv (datafile, header = TRUE, sep = ",", dec = ".")
head(data1,n=20L)
summary(data1)    #what is the highest concentration studied? (can be expressed as both mg/mL and xMIC)

#' example Plot DATA exploration to a PDF file. 
#' 
library(ggplot2)
pdf("DATA.pdf")
plot <- ggplot(data=data1,aes(x=TIME,y=LNDV, colour=as.factor(XMIC))) +
		geom_point()
plot2 <- ggplot(data=data1,aes(x=TIME,y=CAB, colour=as.factor(XMIC))) +
		geom_point()
print(plot)
print(plot2)
dev.off()


# TASK2 END

#' TASK3 start (NONMEN estimation)
#' ==================================================================
#' Estimate (NONMEM) with translated MDL
#+ NONMEN-estimation
cipro <- estimate(cipromdl, "NONMEM", subfolder = sub("^(.+)\\.(.+)$", "\\1", cipromdl))

#' Get population parameters
cipro_params <- getPopulationParameters(cipro, what="estimates")
print("Final estimates:")
print(cipro_params)

#' Xpose diagnostics using NONMEM output
nm.xpdb <- as.xpdb(cipro,datafile)

#' Export results to a PDF file
pdf("GOF_NM.pdf")
 print(basic.gof(nm.xpdb))
 print(ind.plots(nm.xpdb))
dev.off()

# TASK3 END

#' TASK4 start (Strain 202: update inits, VPC and plot)
#' ==================================================================
#' Before running the VPC with PsN we will create a new model file with updated (initial) values in the MDL Parameter Object
#' MLE estimates from previous step can be used 
structuralPar <- getPopulationParameters(cipro, what="estimates",block='structural')$MLE

#' Update the parameter object using the ddmore "updateParObj" function.
myParObj <- getParameterObjects(cipromdl)[[1]]
myParObjUpdated <- updateParObj(myParObj,block="STRUCTURAL",
		item=names(structuralPar),
		with=list(value=structuralPar))

#' Assembling the new MOG. Note that we reuse the data, model and tasks from the previous run.  
myVPCMOG <- createMogObj(dataObj = getDataObjects(cipromdl)[[1]], 
		parObj = myParObjUpdated, 
		mdlObj = getModelObjects(cipromdl)[[1]], 
		taskObj = getTaskPropertiesObjects(cipromdl)[[1]])

#' We can then write the MOG back out to an .mdl file.
ciprofile.VPC <- "cipro202_VPC.mdl"
writeMogObj(myVPCMOG, ciprofile.VPC)

#' Execute VPC with updated model
#+ Strain-202-VPC-execution
vpcFolder <- sub("^(.+)\\.(.+)$", "\\1", ciprofile.VPC)
vpcSO <- VPC.PsN(ciprofile.VPC, command = "vpc", samples=50, seed=12345, vpcOptions ="-stratify_on=CAB -auto_bin=unique", subfolder=vpcFolder, plot=FALSE)

#' Set PDF filename
fname = "VPC_PLOT_202.pdf"

#' Read VPC results
vpc.file <- file.path(vpcFolder, "vpc_results.csv")
vpctab   <- file.path(vpcFolder, "vpctab")
cond.var <- file.path(vpcFolder, "CAB")

#' Open PDF device
pdf(fname, width=11.69, height=8.27)
getwd()

#' Plot VPC through xpose
trellis.par.set(theme=canonical.theme(col=FALSE))
print(xpose.VPC(vpctab=vpctab,vpc.info=vpc.file,PI="both",by="CAB",layout=c(9,1),
          PI.limits=c(0.1,0.9), PI.ci=NULL, PI.real=NULL, PI.area.smooth = T, PI.up.lwd = 1,PI.down.lwd = 1, PI.med.lwd = 1,
          PI.arcol = "grey",col = 1, cex = "0.6", hline =c(2.302609,6.123,13.03),    
          hlcol ="grey60",  hllty =2, hllwd= 1, main=NULL, as.table=T, xlb="Time (h)", ylb="Observed and model predicted bacterial count (log CFU/ml)",
          xlim=c(-3,30), ylim=c(-2,25),
          scales = list(x = list( at=seq(0,24,by=6), 
                        labels=seq(0,24,by=6), relation="same", alternating = c(1), tck=c(1,0), cex="1", tick.number=4),
              y = list(at=c(4.605,9.2103,13.8155,18.4207,23.02585), 
                       labels=c(2,4,6,8,10), relation="same",alternating = c(1), tck=c(1,0), cex="1", tick.number=10))          
))

#' Close PDF device
dev.off()

# TASK4 END

#' TASK5 start (Strain 378: VPC and plot)
#' ==================================================================
#' Set name of modelfile, first you need to you create cipro378_VPC.mdl, see HO instructions.
ciprofile.VPC <- "cipro378_VPC.mdl"

#' Execute VPC with updated model using EC50/MIC correlation
#+ Strain-378-VPC-execution
vpcFolder <- sub("^(.+)\\.(.+)$", "\\1", ciprofile.VPC)
vpcSO <- VPC.PsN(ciprofile.VPC, command = "vpc", samples=50, seed=12345, vpcOptions ="-stratify_on=CAB -auto_bin=unique", subfolder=vpcFolder, plot=FALSE) 


#' Set PDF filename
fname = "VPC_PLOT_378.pdf"

#' Read VPC results
vpc.file <- file.path(vpcFolder, "vpc_results.csv")
vpctab   <- file.path(vpcFolder, "vpctab")
cond.var <- file.path(vpcFolder, "CAB")

#' Open PDF device
pdf(fname, width=11.69, height=8.27)
getwd()

#' Plot VPC through xpose
trellis.par.set(theme=canonical.theme(col=FALSE))
print(xpose.VPC(vpctab=vpctab,vpc.info=vpc.file,PI="both",by="CAB",layout=c(9,1),
				PI.limits=c(0.1,0.9), PI.ci=NULL, PI.real=NULL, PI.area.smooth = T, PI.up.lwd = 1,PI.down.lwd = 1, PI.med.lwd = 1,
				PI.arcol = "grey",col = 1, cex = "0.6", hline =c(2.302609,6.123,13.03),    
				hlcol ="grey60",  hllty =2, hllwd= 1, main=NULL, as.table=T, xlb="Time (h)", ylb="Observed and model predicted bacterial count (log CFU/ml)",
				xlim=c(-3,30), ylim=c(-2,25),
				scales = list(x = list( at=seq(0,24,by=6), 
								labels=seq(0,24,by=6), relation="same", alternating = c(1), tck=c(1,0), cex="1", tick.number=4),
						y = list(at=c(4.605,9.2103,13.8155,18.4207,23.02585), 
								labels=c(2,4,6,8,10), relation="same",alternating = c(1), tck=c(1,0), cex="1", tick.number=10)) 
))

#' Close PDF device
dev.off()

# TASK5 END