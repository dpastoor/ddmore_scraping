;; 1. Based on: Executable_ModelI
;; 2. Description: Morphine PK across paediatric age-range
;; x1. Author: user
$PROB Morphine PK across paediatric age-range
$INPUT 
ID 
TIME ;in min
AMT  ;in microgram (ug)
RATE ;in ug/min
DV   ;natural logarithm of concentration
MDV 
CONC ;ug/L
CMT  ;1=morphine central plasma
BW   ;bodyweight in kg 
POP  ;1 = 0 - 3 years, 2 = 6 - 15 years, 3= 18 - 36 years (see publication)

$DATA Simulated_DataModel1.csv IGNORE=I
$SUBROUTINES ADVAN5
;-------------------------------------------------------------------------------
$MODEL
NCOMPARTMENTS=2  
COMP (CENTRAL, DEFDOSE)  ;MORPHINE CENTRAL
COMP =(2)                ;PERPHERAL COM OF MORHINE
;-------------------------------------------------------------------------------
$PK
KDEC    = THETA(1)                      ; DECREASE OF EXPONENT FOR CLM1
KMAX    = THETA(2)+KDEC                 ; MAXIMUM EXPONENT OF CLM1
KHAL    = THETA(3)                      ; K50 OF CLM1
GAMMA   = THETA(4)                      ; GAMMA OF CLM1
KBDE    = KMAX-KDEC*(BW**GAMMA)/(KHAL**GAMMA+BW**GAMMA)

;-------------------------------------------------------------------------------
TVCL    = THETA(5)*(BW/70)**KBDE        ; POPULATION CLEARANCE OF MORPHINE
CL      = TVCL*EXP(ETA(1))              ; INDIVIDUAL ...
TVQ2    = THETA(6)*(BW/70)              ; POPULATION INTERCOMPARTMENTAL CLEARANCE OF MORHPINE
IF (POP.EQ.2) TVQ2= THETA(10) 
Q2      = TVQ2*EXP(ETA(2))              ; INDIVIDUAL ...

;-------------------------------------------------------------------------------
TVV1    = THETA(7)*(BW/70)              ; POPULATION VOLUME OF MORPHINE CENTRAL COMPARTMENT
IF (POP.EQ.2) TVV1=THETA(11)*(BW/70);
V1      = TVV1 * EXP(ETA(3))            ; INDIVIDUAL ...

TVV2    = THETA(8)*(BW/70)              ; POPULATION VOLUME OF MORPHINE PEREPHERAL COMPARTMENT
V2      = TVV2*EXP(ETA(4))              ; INDIVIDUAL ...

;--------------------------------------------------------------------------------------------
ET1  = ETA(1)
ET2  = ETA(2)
ET3  = ETA(3)
ET4  = ETA(4)

;-------------------------------------------------------------------------------
S1=V1
;-------------------------------------------------------------------------------
K10  = CL/V1
K12  = Q2/V1
K21  = Q2/V2 

F1      = 1
IF (POP.EQ.3)  F1 = 0.88

;-------------------------------------------------------------------------------
$ERROR
IPRED=LOG(0.000001)
IF (F.GT.0) IPRED = LOG(F) 
W   =  THETA(9)
IRES  = IPRED-DV
IWRES = IRES/W
Y = IPRED + ERR(1)*W    

;-------------------------------------------------------------------------------
$THETA
(0.1,  0.594,   )        ;KDEC (TH1)
(0.2,  0.872,   )        ;KMAX-KDEC OR MINMUM EXP (TH2)
(0.05, 4.01, 20 )        ;KHAL (TH3)
(1,    4.62,    )        ;GAMMA (TH4)
(0.001,1.62,    )        ;CL (TH5)
(0.01, 1.90,    )        ;Q2 (TH6)
(0.1, 81.2,     )        ;V1 (TH7)
(0.1,128,       )        ;V2 (TH8)
(0,    0.432,   )        ;ERR1 (TH9)
(0.1,  0.500,   )        ;Q2 ADO (TH10)
(0.1, 46.0,     )        ;V1 ADO (TH11)
;-------------------------------------------------------------------------------
$OMEGA
0.159    ;CL
0 FIX    ;Q2
0.253    ;V1
0 FIX    ;V2
;-------------------------------------------------------------------------------
$SIGMA
1  FIX   ; ERR1
;-------------------------------------------------------------------------------

$EST NOABORT SIGDIG=3 PRINT=15 MAXEVAL=9999 METHOD=1 INTERACTION POSTHOC
$COV COMP PRINT=E
$TABLE ID TIME AMT RATE DV PRED IWRES CWRES IPRED MDV CONC CMT BW POP KBDE TVCL CL Q2 TVV1 V1 V2 ET1 ET2 ET3 ET4
NOPRINT ONEHEADER NOAPPEND FILE=sdtabSimulatedDataModel1
