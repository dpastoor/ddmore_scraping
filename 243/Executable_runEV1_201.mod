


$PROBLEM Gompertz hazard model. runEV1_102.    

$INPUT ID AGE NEUT PRE_TRE MAX_LEG AUC CMAX DBLOCK TIME DV DVID EVID 

$DATA Simulated_event_data.csv
ACCEPT=(DVID.EQ.1,DVID.EQ.0)                                                        ; Select which DVID to accept

$SUBROUTINE ADVAN13 TRANS1 TOL=9

$MODEL NCOMPARTMENTS=1
       COMP=(CPT1)      	; TTE model			



$PK 

CENSORING = 1

IF (NEWIND.LE.1) SURVZ  = 1   		                ; Survival(0)=1 [we use this variable for storing survival function at start of observation interval for interval censored data]




;;;;;;;;;;;;;;;;;;;;; Start of TTE model ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Lam          = THETA(1)*EXP(ETA(1)) 		; NONMEM needs an eta, but it will be fixed to zero	
	
			
LamC   = Lam/1000				; rescale the parameters so NONMEM does not have to search tiny numbers


VAL          = EXP((THETA(2)/10000)*(NEUT-4133))*EXP((THETA(3)/100)*(AGE-55))




$DES
DADT(1) = VAL*LamC



$ERROR

CHAZ = A(1)				; cumulative hazard
SURV = EXP(-CHAZ)  			; probability of surviving to or beyond current time


HAZNOW = VAL*LamC


                                                                  


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
(0,4.7,50000)
(-1000,-0.1,1000)
(-1000,3,1000)


$OMEGA
0 FIX 
	

$EST MAXEVAL=9999 METHOD=COND LAPLACE NUMERICAL LIKE SLOW NOABORT NOTHETABOUNDTEST NSIG=4 SIGL=9 PRINT=5 MSFO=runEV1_201.txt

$COV PRINT=E 

$TABLE ID AGE NEUT PRE_TRE MAX_LEG AUC CMAX DBLOCK TIME DV DVID EVID
SURV CHAZ HAZNOW FORMAT=s1PE18.9 NOPRINT ONEHEADER FILE=runEV1_201.tab

