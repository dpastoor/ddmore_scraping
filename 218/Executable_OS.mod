; TTE model for OS. Weibull baseline hazard. Covariates are SLD at baseline, TSR(t) for t<= week 12 and TSR(week12) for t> week 12, appearance of new lesion (dichotomous time-varying covariate) and ECOG status at enrolment
; Parameter of the SLD(t) model were estimated previously (IPP approach)

$SIZES LIM6=2000
$PROBLEM OS-IPP

$INPUT CCOM,
ID,
TIME, ;day
DV, ;0 if censored, 1 if death 
CB, ;exposure to carboplatin (per-cycle average AUC)
G, ;exposure to gemcitabine (per-cycle average AUC)
EVID,
FLG, ;FLG=9 for OS data, FLG=1 for covariates-related entries
CMT, ;CMT=1 for SLD, CMT=2 for hazard of death
NWLS, ;dichotomous time-varying covariate indicating the presence of new lesions with respect to enrolment
KG, ;subjec specific tumour growth rate constant (estimated previously)
KD0, ;subjec specific carboplatin related tumour death rate constant (estimated previously)
KD1, ;subjec specific gemcitabine-related tumour death rate constant (estimated previously)
IBASE, ;subjec specific SLD at baseline (estimated previously)
SLD0, ;measured SLD at enrolment
ECOG, ;measured ECOG status at enrolment
OSCENS ;OSCENS=0 if patient died, OSCENS=1 if censoring occured

$DATA Simulated_OS.csv IGNORE=C 

$SUBROUTINE ADVAN6 TOL=3

$MODEL
COMP = (TUMOUR) 
COMP = (HDEATH)

$PK
; --- OS param ---
LAM = THETA(1)*EXP(ETA(1)) ;scale parameter
SHP = THETA(2) ;shape parameter
BSLD0 = THETA(3) ;parameter for SLD at enrolment
BTSR = THETA(4) ;parameter for TSR(t) (tumour size ratio)
BNWLS = THETA(5) ;parameter for NewLesion(t)
BECOG = THETA(6) ;parameter for ECOG at enrolment


; ==== SLD baseline ====
A_0(1) = IBASE*1000
MMBAS = IBASE*1000

; --- Bacward interpolation of covariates
; exposure to drug
IF(NEWIND.NE.2) OCB=CB
IF(NEWIND.NE.2) OG=G
E0 = OCB ;ng/dL*day/n days in cycle
E1 = OG ;mol/10^6 cells*day/n days in cycle ;sum Parent and active metabolite
OCB=CB
OG=G
; Appearance of new lesion
IF(NEWIND.NE.2) ONWLS=NWLS
INWLS=ONWLS
ONWLS=NWLS

; --- Time constant covariates
; Normalised SLD at enrolment
TVSLD0 = 70 ;average SLD at enrolment
NSLD0 = SLD0/TVSLD0 ;normalised SLD at enrolment
; ECOG at enrolment
IECOG=ECOG

$DES
DADT(1) = KG/1000 * A(1) - (KD0/1000 * E0 + KD1/100 * E1) * A(1) 
TUM = A(1)
TSR = (TUM-MMBAS)/MMBAS
; --- OS ----
IF(T.EQ.0) THEN
	WTS=0
	TM12=0
ENDIF
IF(T.LE.84) THEN
WTS = TSR
TM12 = WTS
ELSE
WTS=WTS
ENDIF

DEL = 1E-6
DADT(2) = LAM*SHP*(LAM*(T+DEL))**(SHP-1) *EXP(BSLD0*NSLD0+BTSR*WTS+BNWLS*INWLS+BECOG*IECOG)

$ERROR
DELX = 1E-6
XTUM=A(1)
XTSR = (XTUM-MMBAS)/MMBAS
IF(TIME.EQ.0) THEN
	XWTS=0
	XTM12=0
ENDIF
IF(TIME.LE.84) THEN
	XWTS = XTSR
	XTM12=XWTS
ELSE
	XWTS=XWTS
ENDIF

;--- Death hazard ---
CHZ = A(2)
SUR=EXP(-CHZ)

HAZN = LAM*SHP*(LAM*(TIME+DELX))**(SHP-1)*EXP(BSLD0*NSLD0+BTSR*XWTS+BNWLS*INWLS+BECOG*IECOG)

IF (FLG.EQ.9.AND.EVID.EQ.0.AND.OSCENS.EQ.1) THEN
IPRED=SUR ;probability of survival (censored event)
Y = IPRED  ; Y is probability for TTE data
ENDIF
IF (FLG.EQ.9.AND.EVID.EQ.0.AND.OSCENS.EQ.0) THEN
IPRED=SUR*HAZN ;probability of event (death) at time=TIME
Y = IPRED  ; Y is probability for TTE data
ENDIF

;Parameters
$THETA
; --- 0S ---
(0, 0.001) ;1 LAM
(0, 2) ;2 SHP

(0, 0.1) ;3 BSLD0
(0, 0.1) ;4 BTSR
(0, 0.1) ;5 BNWLS
(0, 0.1) ;6 BECOG

$OMEGA 
 0 FIX ; KG ;placeholder

$EST MAXEVAL=9000 PRINT=10 METH=COND LAPLACE LIKE NOABORT NSIG=3 SIGL=9
$COV SLOW

$TABLE ID TIME DV IPRED EVID FLG XWTS 
FILE=sdtab_OS NOPRINT ONEHEADER


