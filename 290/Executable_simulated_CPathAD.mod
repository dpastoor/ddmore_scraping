;; 1. Based on: run47
;; 2. Description: run14.mod CAMD ADAS-COG TOTAL
;; x1. Author: Daniela Conrado
;;Tue Apr  8 11:48:30 EDT 2014
$PROBLEM run48.mod CAMD ADAS-COG TOTAL

$INPUT
ID	
STUDY	
TIME	
WEEK	
ADAS	
DROPOUT	
MMSE	
AGE	
SEX	
DXDATE	
START	
END	
QSDTC	
APOE4	
COMED	
COMEDEND	
COMEDSTART	
COMEDENDIND	
COMEDSTARTIND	
PRIMCOMED	
PRIMCOMEDS	
ADASBL	
MDV	
ADASMOD	
ADASTRANS=DV	
COMED2	
APOE4C

;ID: patient identification
;STUDY: study identification
;VISIT: study visit
;WEEK: study week
;TIME: time (days)
;ADAS: ADAS-cog total score (0-70)
;DROPOUT: 0=patient has never dropped out; 1=patient has dropped out at some point
;It does not directly carry time information, meaning that each subject has been assigned to one dropout code
;Do not know when exactly the patient dropped out but can infer based on the last time of the non-missing DV
;COVARIATES AT MMSE BASELINE:
;MMSE: BASELINE mini-mental state evaluation (normally taken at screening)
;AGE: age (years)
;SEX: 1=female, 2=male
;DXDATE: year of diagnosis
;RACE: 1="AMERICAN INDIAN OR ALASKA NATIVE", 2="ASIAN", 3="BLACK OR AFRICAN AMERICAN", 4/5="MULTIPLE"/"MUTIPLE", 
;6="NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", 7="OTHER", 8="UNKNOWN", 9="WHITE"
;FAMILY: family history of AD 1=no, 2=unknown, 3=yes
;START: patient start year
;END: patient end year
;QSDTC: CDISC standard term for Date/Time of Finding
;APOE: 0=non-carrrier, 1=heterozygous, 2=homozygous, 3=unknown
;MDV: missing dependent variable
;ADASMOD: (ADAS * (nid-1) + 0.5)/nid, where nid is the total number of individuals
;Based on: Smithson & Verkuilen. Psychological Methods, 11: 54-71, 2006 [please refer to page 55]
;ADASTRANS: (ADASMOD - 0)/(70 - 0)
;ADASBLT2TRANS: adas baseline including the ones with baseline indicator; for missing bl, adas at time 0 or 1 was used in this order
;Based on: Smithson & Verkuilen. Psychological Methods, 11: 54-71, 2006 [please refer to page 54]
;INITIAL ESTIMATES:unless stated differently, it was based on lme model for ADASTRANS in R
;MOTIVATION FOR USING BETA-REGRESSION: for bounded data, the expectation must be nonlinear due to the ceiling/floor effects,
;and the error distribution must be heteroscedastic since the variance must approach zero as their mean approaches
;either boundary score (Reference: Xu et al. J PK PD, 40: 537-544, 2013)

$DATA ./r_code_simulation_ddmore/simulated_dataset.csv IGNORE=@ IGNORE=(MDV.EQ.1,ID.EQ.528);

$OMEGA  BLOCK(2)
 0.156  ;   1 BSV_BL
 0.0224  ; CORR BSV_BL:BSV_SL
 0.0424 ; 2 BSV_SL

$OMEGA  BLOCK(2)
0.0084 ;   1 ISV_BL
0.0014 ;   2 CORR ISV_BL:ISV_SL
0.0004 ;   3 ISV_SL

;##########################################################################################################################################

$PRED

;;; BLAPOE4C-DEFINITION START
BLAPOE4C = ( 1 + THETA(9)*(APOE4C - 0.72))
;;; BLAPOE4C-DEFINITION END


;;; SLCOMED2-DEFINITION START
IF(COMED2.EQ.1) SLCOMED2 = 1  ; Most common
IF(COMED2.EQ.0) SLCOMED2 = 1+THETA(8)
;;; SLCOMED2-DEFINITION END


;;; SLAPOE4C-DEFINITION START
SLAPOE4C = ( 1 + THETA(7)*(APOE4C - 0.72))
;;; SLAPOE4C-DEFINITION END


;;; SLAGE-DEFINITION START
SLAGE = ( 1 + THETA(6)*(AGE - 75.00))
;;; SLAGE-DEFINITION END

;;; SL-RELATION START
SLCOV=SLAGE*SLAPOE4C*SLCOMED2
;;; SL-RELATION END


;;; BLSEX-DEFINITION START
IF(SEX.EQ.1) BLSEX = 1  ; Most common
IF(SEX.EQ.2) BLSEX = THETA(5)
;;; BLSEX-DEFINITION END

;;; BL-RELATION START
BLCOV=BLSEX*BLAPOE4C
;;; BL-RELATION END

YTIME	 = TIME/365.25									;TIME in years to avoid such small slope
;############################### DEFINING THETAS ###########################################################################################

TVBL     = THETA(1)				 					    ;POPULATION BASELINE ADAS-Cog

TVBL = BLCOV*TVBL
TVSL     = THETA(2)                					    ;POPULATION SLOPE: disease progression rate constant

TVSL = SLCOV*TVSL


STEWSV = 1
STEBSV = 1
IF(STUDY.EQ.1131) THEN
  STEWSV = EXP(THETA(10))
  STEBSV = EXP(THETA(11))
ENDIF

BL       = TVBL*EXP(ETA(1)*STEBSV + ETA(3))        			    ;EXPONENTIAL ERROR: to prevent individual baseline from becoming negative
SL       = TVSL + ETA(2) + ETA(4)		     			;ADDITIVE ERROR: slope can be positive or negative on individual basis

TVSHAPE  = THETA(3)				 					    ;POPULATION SHAPE factor of the Richards logistic growth model
TAU      = THETA(4)*STEWSV				 					    ;TAU is the precision parameter of the beta distribution (=ALPHA+BETA)
														;VARIANCE of the beta distribution is MUR*(1-MUR)/(TAU+1), 
														;meaning that the variance increases as TAU increase; although, it is not
														;the sole determinante of dispersion, since dispersion also depends on MUR

;############################### STRUCTURAL MODEL TO DESCRIBE THE RATE OF DISEASE PROGRESSION ##############################################
														;Using Richards logistic growth model: three-parameter logistic model
														;Tsoularis & Wallace. Mathematical Biosciences, 179:21-55, 2002 [please refer to page 33]

;GG=1.0E-30
;IF(BL<GG) THEN
;   BL=GG
;   PRINT 1,ID
;ENDIF
DEN1     = BL**TVSHAPE
DEN2     = (70**TVSHAPE) - (BL**TVSHAPE)
DEN3     = EXP(-TVSHAPE*SL*YTIME)
DENN=DEN1 + DEN2*DEN3
;IF(DENN<GG) THEN
;  DENN=GG
;  PRINT 2, ID
;ENDIF
MUR      = BL/((DENN)**(1/TVSHAPE)) 	     			;INDIVIDUAL ADAS-Cog response as a function of time


;############################### DEFINING IPRED/PRED ########################################################################################

F        = MUR					      					;INDIVIDUAL (IPRED)
IPRED    = F

;############################## RATE OF CHANGE IN ADAS-Cog AS A FUNCTION OF TIME ############################################################

;RATE     = SL*MUR*(1-((MUR/70)**TVSHAPE)) 				;as differential equation: RATE=d(ADAS-Cog)/dt
;TVRATE   = TVSL*TVMUR*(1-((TVMUR/70)**TVSHAPE))		    ;population value

;############################## BETA REGRESSION #############################################################################################
														;Reference: Smithson & Verkuilen. Psychological Methods, 11: 54-71, 2006 [page 58]
														;MUR is the location parameter of the beta distrution
														;TAU is the precision parameter of the beta distribution
														;MUR and TAU have no restriction on each other, so that they can be modeled separately
ALPHA    = MUR*TAU										;BETA DISTRIBUTION shape parameter pulling density toward 0 (alpha>0)
BETA     = (1-MUR)*TAU								    ;BETA DISTRIBUTION shape parameter pulling density toward 1 (beta >0)

X1       = ALPHA+BETA									;It is really TAU
X2       = ALPHA
X3       = BETA

;NEMES APPROXIMATION OF THE LN(GAMMA) DISTRIBUTION
GG1=1.0E-150
;IF(X1<GG1) THEN
;   X1=GG1
;   PRINT 3, ID
;ENDIF
;IF(X2<GG1) THEN
;  X2=GG1
;  PRINT 4,ID
;ENDIF
;IF(X3<GG1) THEN
;  X3=GG1
;  PRINT 5,ID
;ENDIF
LGAMMAX1 = X1*(LOG(X1)-1) + 0.5*(LOG(2*3.1415)-LOG(X1)) + (5/4)*X1*LOG(1+(1/(15*X1**2)))
LGAMMAX2 = X2*(LOG(X2)-1) + 0.5*(LOG(2*3.1415)-LOG(X2)) + (5/4)*X2*LOG(1+(1/(15*X2**2)))
LGAMMAX3 = X3*(LOG(X3)-1) + 0.5*(LOG(2*3.1415)-LOG(X3)) + (5/4)*X3*LOG(1+(1/(15*X3**2)))

;LOG-LIKELIHOOD OF THE BETA DISTRIBUTION: LN of the probability density function of the beta distribution
WDV=DV
IF(WDV<0.001) THEN
   WDV=0.001
;   PRINT 6,ID
ENDIF
IF(WDV>0.999) THEN
   WDV=0.999
;   PRINT 7,ID
ENDIF
LLBETA   = LGAMMAX1 - LGAMMAX2 - LGAMMAX3 + (ALPHA-1)*LOG(WDV) + (BETA-1)*LOG(1-WDV)

;############################### DEFINING RESIDUALS AND Y ###################################################################################

;PEARSON RESIDUALS
;SOR      = (WDV-IPRED)/SQRT(IPRED*(1-IPRED)/(1+TAU))    ;Standardized ordinary residuals

LLB2=-2.0*LLBETA
;IF(LLB2.GT.300.0) THEN
;   LLB2=300.0
;   PRINT 8,ID
;ENDIF 
;IF(LLB2.LT.-300.0) THEN
;   LLB2=-300.0
;   PRINT 9,ID
;ENDIF
Y        = LLB2
;############################## INITIAL ESTIMATES ###########################################################################################

$THETA
(0, 22.3, 100) 				; 1 TVBL
(-100, 0.151, 50) 			; 2 TVSL
(0, 6.98, 100) 				; 3 TVSHAPE (initial estimate from bapineuzumab model)
(-300, 87.8, 300) 		        ; 4 TAU (initial estimate from bapineuzumab model)
(0,1)                            	; 5 BLSEX1
(0.1)                   		; 6 SLAGE1
(0.2)    		                ; 7 SLAPOE4C1
(-0.5, -0.275, 0.1)                   		; 8 SLCOMED21
(0.1)  	                                ; 9 BLAPOE4C1
(0.1)                                   ; 10 STE WSV EFFECT
(0.1)                                   ; 11 STE BSV EFFECT


$LEVEL
STUDY=(3[1],4[2])

$ESTIMATION MAXEVAL=9999 PRINT=1 METHOD=BAYES -2LL LAPLACIAN NOHABORT FNLETA=0 NSIG=2 SIGL=7 SLOW MCETA=10 NONINFETA=1 FILE=run48.ext
;$ESTIMATION MAXEVAL=9999 PRINT=1 METHOD=COND -2LL LAPLACIAN NOHABORT SIGL=10 FILE=run48.ext

$COV MATRIX=R UNCONDITIONAL SIGL=9 SLOW PRINT=E

$TABLE STUDY ID TIME ADAS WDV MUR NOPRINT ONEHEADER FILE=sdtab48.txt
;$TABLE STUDY ID TIME ADAS RATE TVRATE WDV MUR TVMUR SOR NOPRINT ONEHEADER FILE=sdtab48.txt
;$TABLE STUDY ID BL SL TVSHAPE ALPHA BETA TAU ETA(1) ETA(2) ETA(3) ETA(4) NOPRINT ONEHEADER FIRSTONLY FILE=patab48.txt
$TABLE STUDY ID BL SL TVSHAPE ALPHA BETA TAU ETA(1) ETA(2) NOPRINT ONEHEADER FIRSTONLY FILE=patab48.txt
$TABLE STUDY ID AGE SEX START APOE4 APOE4C COMED2 MMSE BL ADASBL NOPRINT ONEHEADER FILE=cotab48.txt