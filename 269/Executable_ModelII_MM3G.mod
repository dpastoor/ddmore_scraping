$PROB Morphine and M3G PK across paediatric age-range
$INPUT 
ID 
TIME ;in min 
AMT  ;in microgram (ug)
RATE ;ug/min
DV   ;natural logarithm of concentration
MDV 
CONC ;in ug/L
CMT  ;1=morphine central plasma, 2=M3G in plasma
BW   ;bodyweight in kg
POP  ;value = 2 --> young children - adolscents, value !=2 --> newborns and young children and adults (see publication)

$DATA Simulated_DataModel2.csv IGNORE=I 
$SUBROUTINES ADVAN5
;-------------------------------------------------------------------------------
$MODEL
NCOMPARTMENTS=3  
COMP (CENTRAL, DEFDOSE) ;MORPHINE CENTRAL
COMP=(2)                ;M3G
COMP=(3)                ;PERPHERAL COM OF MORPHINE
;-------------------------------------------------------------------------------
$PK
KDEC1   = THETA(1)                      ; DECREASE OF EXPONENT FOR CLM1
KMAX1   = THETA(2)+KDEC1                ; MAXIMUM EXPONENT OF CLM1
KHAL1   = THETA(3)                      ; K50 OF CLM1
GAMMA1  = THETA(4)                      ; GAMMA OF CLM1
KBDE1   = KMAX1-KDEC1*(BW**GAMMA1)/(KHAL1**GAMMA1+BW**GAMMA1)

KDEC2   = THETA(5)                      ; DECREASE OF EXPONENT FOR CLE2
KMAX2   = THETA(6)+KDEC2                ; MAXIMUM EXPONENT OF CLE2
KHAL2   = THETA(7)                      ; K50 OF CLE2
GAMMA2  = THETA(8)                      ; GAMMA OF CLE2
KBDE2   = KMAX2-KDEC2*(BW**GAMMA2)/(KHAL2**GAMMA2+BW**GAMMA2)

;-------------------------------------------------------------------------------
TVCLM1  = THETA(9)*(BW/70)**KBDE1       ; POPULATION METABOLISM OF MORPHINE TO M3G
CLM1    = TVCLM1*EXP(ETA(1))            ; INDIVIDUAL ...
TVCLE1  = THETA(10)*(BW/70)             ; POPULATION EXCRETION OF MORPHINE + METABOLISM TO M6G
CLE1    = TVCLE1*EXP(ETA(2))            ; INDIVIDUAL ...
TVCLE2  = THETA(11)*(BW/70)**KBDE2      ; POPULATION EXCRETION OF M3G
CLE2    = TVCLE2*EXP(ETA(3))            ; INDIVIDUAL ...
;-------------------------------------------------------------------------------
TVV1    = THETA(12)*(BW/70)             ; POPULATION VOLUME OF MORPHINE CENTRAL COMPARTMENT
V1      = TVV1 * EXP(ETA(4))            ; INDIVIDUAL ...
TVV2    = THETA(13)*(BW/70)**THETA(18)  ; POPULATION VOLUME OF M3G
V2      = TVV2 * EXP(ETA(5))            ; INDIVIDUAL ...
TVQ2    = THETA(14)*(BW/70)             ; POPULATION INTERCOMPARTMENTAL CLEARANCE OF MORHPINE
Q2      = TVQ2*EXP(ETA(6))              ; INDIVIDUAL ...
TVV3    = THETA(15)*(BW/70)             ; POPULATION VOLUME OF MORPHINE PEREPHERAL COMPARTMENT
V3      = TVV3*EXP(ETA(7))              ; INDIVIDUAL ...

F1      = 1
IF (POP.EQ.3) F1 = 0.88
;-------------------------------------------------------------------------------
ET1  = ETA(1)
ET2  = ETA(2)
ET3  = ETA(3)
ET4  = ETA(4)
ET5  = ETA(5)
ET6  = ETA(6)
ET7  = ETA(7)
;-------------------------------------------------------------------------------
S1=V1
S2=V2
;-------------------------------------------------------------------------------
K10  = CLE1/V1
K12  = CLM1/V1
K13  = Q2/V1
K20  = CLE2/V2
K31  = Q2/V3

;-------------------------------------------------------------------------------
$ERROR
COM1=0
IF (CMT.EQ.1) COM1=1
COM2=0
IF (CMT.EQ.2) COM2=1

IPRED=LOG(0.000001)
IF (F.GT.0) IPRED = LOG(F)
W1   =  THETA(16)          ;ERR Morphine
W2   =  THETA(17)          ;ERR M3G

IRES  = IPRED-DV
IWRES = IRES/(COM1*W1+COM2*W2)

Y1 = IPRED + ERR(1)*W1     ; MORPHINE
Y2 = IPRED + ERR(2)*W2     ; M3G

Y=COM1*Y1+COM2*Y2   
;-------------------------------------------------------------------------------
$THETA
(0.1,  0.665   )        ;KDEC1
(0.4,  0.890   )        ;KMAX1-KDEC1 OR MINMUM EXP1 (TH2)
(0.05, 3.89,20 )        ;KHAL1
(1,    3.61    )        ;GAMMA1
;--------------------------------------------------------
(0.1,  0.448   )        ;KDEC2
(0.4,  0.610   )        ;KMAX2-KDEC2 OR MINMUM EXP2 (TH15)
(0.05, 4.87, 20)        ;KHAL2
(1,    6.84    )        ;GAMMA2
;---------------------------------------------------------
(0.001, 1.67   )        ;CLM1
(0.001, 0.0572 )        ;CLE1
(0.001, 0.225  )        ;CLE2
(0.01 , 29.3   )        ;V1
20 FIX                   ;V2
(0.001, 4.20   )        ;Q2
(0.01,155      )        ;V4 
(0,     0.447  )        ;ERR MORP
(0,     0.371  )        ;ERR M3G
(0.1,   0.711  )        ;EXP V2
;-------------------------------------------------------------------------------
$OMEGA
0.202     ;CLM1
0.0674    ;CLE1
0.191     ;CLE2
0.508     ;V1
0.367     ;V2
0 FIX     ;Q2
0.308     ;V3
;-------------------------------------------------------------------------------
$SIGMA
1 FIX
1 FIX
;-------------------------------------------------------------------------------
$EST NOABORT SIGDIG=3 PRINT=15 MAXEVAL=9999 METHOD=1 INTERACTION POSTHOC
$COV COMP PRINT=E
$TABLE ID TIME AMT RATE DV PRED CWRES IWRES IPRED MDV CONC CMT BW POP TVCLM1 KBDE1 KBDE2 CLM1 TVCLE1 CLE1 TVCLE2 CLE2 Q2 TVV1 V1 TVV2 V2 TVV3 V3 ET1 ET2 ET3 ET4 ET5 ET6 ET7
NOPRINT ONEHEADER NOAPPEND FILE=sdtabSimulatedDataModel2

