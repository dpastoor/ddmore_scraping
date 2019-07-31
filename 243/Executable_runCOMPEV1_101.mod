

$PROBLEM Log_normal model. runCOMPEV1_101.    

$INPUT ID AGE NEUT PRE_TRE MAX_LEG AUC CMAX DBLOCK TIME DV DVID EVID
AGE_CAT NEUT_CAT PRE_TR_C MAX_CAT AUC_CAT CMAX_CAT

$DATA Simulated_event_data.csv
ACCEPT=(DVID.EQ.3,DVID.EQ.0)                                                        ; Select which DVID to accept

$SUBROUTINE ADVAN13 TRANS1 TOL=9

$MODEL NCOMPARTMENTS=1
       COMP=(CPT1)      	; TTE model			




$PK 

CENSORING = 2

IF (NEWIND.LE.1) SURVZ  = 1   		                ; Survival(0)=1 [we use this variable for storing survival function at start of observation interval for interval censored data]


;;;;;;;;;;;;;;;;;;;;; Start of TTE model ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

lambda      = THETA(1)*exp(ETA(1))
alpha	    = THETA(2)


pi  = 3.1415927
fac = (2*pi)**0.5

DEL = 1E-8


VAL          = EXP((THETA(3)/1000)*(AGE-55))


$DES 

num = log((T+DEL)/alpha)

pdf = (1.0/(fac*lambda*(T+DEL)))*exp(-(num**2)/(2*lambda*lambda))
DADT(1)  = VAL*pdf/(1-phi(num/lambda))               			            



$ERROR

CHAZ = A(1)				; cumulative hazard
SURV = EXP(-CHAZ)  			; probability of surviving to or beyond current time

num1 = log((TIME+DEL)/alpha)

pdf1 = (1.0/(fac*lambda*(TIME+DEL)))*exp(-(num1**2)/(2*lambda*lambda))
HAZNOW = VAL*pdf1/(1-phi(num1/lambda))                    			            
             			   



;;;;;;;;;; Right censoring only ;;;;;;;;;;;;;;;;;

CS=0
IF (CENSORING.EQ.1.AND.DV.EQ.-1.AND.TIME.GT.0) CS=1                                  ; CS=1 for right censored events
IF (CENSORING.EQ.1.AND.TIME.GT.0) Y = (1 - CS)*HAZNOW*SURV + CS*SURV                 ; for right censored events the likelihood of event at this time or greater is the survival function
					                                             ; for non-censored events the likelihood of event at this time is (hazard function)*(survival function)
;;;;;;;;;; End ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;; Right + interval censoring ;;;;;;;;;;;;;;;;;;;;;;

IF (CENSORING.EQ.2.AND.DV.EQ.-1) Y=SURV

IF (CENSORING.EQ.2.AND.DV.EQ.2) SURVZ=SURV
	 
IF (CENSORING.EQ.2.AND.DV.NE.2) SURVZ=SURVZ

IF (CENSORING.EQ.2.AND.DV.EQ.3) Y=SURVZ-SURV

;;;;;;;;;; End ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
				

$THETA
(0,1 )
(0,1 )
(-10000,1,10000)


$OMEGA
0 FIX 
	

$EST MAXEVAL=9999 METHOD=COND LAPLACE NUMERICAL LIKE SLOW NOABORT NOTHETABOUNDTEST NSIG=3 SIGL=9 PRINT=5 MSFO=runCOMPEV1_101.txt

$COV PRINT=E 

$TABLE ID AGE NEUT PRE_TRE MAX_LEG AUC CMAX DBLOCK TIME DV DVID EVID
AGE_CAT NEUT_CAT PRE_TR_C MAX_CAT AUC_CAT CMAX_CAT
SURV CHAZ HAZNOW FORMAT=s1PE18.9 NOPRINT ONEHEADER FILE=runCOMPEV1_101.tab


