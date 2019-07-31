$PROB morphine PK in children < 3 years 
$INPUT ID TIME AMT RATE DV MDV CMT BWS PNA SGA  
$DATA Simulated_PaediatricMorphinePK.csv IGNORE=@ 

; TIME in minutes
; AMT in micrograms
; DV log-transformed (natural logarithm) in micrograms per liter
; CMT 1 = morphine, 2 = M3G, 3 = M6G
; BWS bodyweight in grams at time of study
; PNA postnatal age in days
; SGA small for gestational age 1 = yes, 0 = no
; NKOD for extra additive error on 24 hr post-infusion sample in study of Nancy

$SUBROUTINE ADVAN5
$MODEL

COMP (CENTRAL, DEFDOSE)
COMP=(2)                ;M3G
COMP=(3)                ;M6G
COMP (PERIPH)           ;Peripheral cmt for parent

$PK
CL2=(THETA(1)*(BWS*0.001)**THETA(4))*EXP(ETA(1))                  ;parent to M3G (CMT2)
IF (PNA.GT.10) CL2=(THETA(10)*(BWS*0.001)**THETA(4))*EXP(ETA(1)) 
CL2N=(THETA(1)*(BWS*0.001)**THETA(4))
CL2P=(THETA(10)*(BWS*0.001)**THETA(4))

V1 =(THETA(2)*(BWS*0.001)**THETA(5))*EXP(ETA(2))                  ;V1 parent (CMT1)
V1P=(THETA(2)*(BWS*0.001)**THETA(5))
Q  = THETA(3)                                                     ;parent CMT1 to CMT4

CL3=(THETA(6)*(BWS*0.001)**THETA(4))                              ;parent to M6G (CMT3)
IF (PNA.GT.10) CL3=(THETA(11)*(BWS*0.001)**THETA(4))
CL3N=(THETA(6)*(BWS*0.001)**THETA(4))
CL3P=(THETA(11)*(BWS*0.001)**THETA(4))

CL4= (THETA(7)*(BWS*0.001)**THETA(4))*EXP(ETA(4))                 ;excretion M3G 
CL5= (THETA(9)*(BWS*0.001)**THETA(4))*EXP(ETA(3))                 ;excretion M6G 

V2 = V1*THETA(8)                                                  ;V M3G
V3 = V2                                                           ;V M6G
V4 = V1                                                           ;V parent peripheral

K12=CL2/V1
K13=CL3/V1
K20=CL4/V2
K30=CL5/V3
K14=Q/V1
K41=Q/V4

S1=V1
S2=V2
S3=V3

ET1=ETA(1)
ET2=ETA(2)
ET3=ETA(3)
ET4=ETA(4)

$ERROR
COM1=0
IF (CMT.EQ.1) COM1=1
COM2=0
IF (CMT.EQ.2) COM2=1
COM3=0
IF (CMT.EQ.3) COM3=1

IPRE=LOG(0.000001)
IF (F.GT.0) IPRE = LOG(F)
Y1=IPRE+ERR(1)
Y2=IPRE+ERR(2)
Y3=IPRE+ERR(3)
Y=COM1*Y1+COM2*Y2+COM3*Y3
LIPR=F
IWRE=0
IF (IPRE.GT.0) IWRE=(DV-IPRE)/IPRE

$THETA
(0, 0.00309) ; TH1 CL2 <10days
(0, 1.99) ; TH2 V1
(0, 0.0289) ; TH3 Q
(0, 1.44) ; TH4 exponent BW on CL
(1) FIX ; TH5 exponent BW on V
(0, 0.00041) ; TH6 CL3 <10 days
(0, 0.00219) ; TH7 CL4
(0, 0.119) ; TH8 V3 and V4 (fraction of V1)
(0, 0.00111) ; TH9 CL5
(0, 0.00825) ; TH10 CL2 >10days
(0, 0.0007) ; TH11 CL3 >10days

$OMEGA
 0.104 ; OM1 CL2
 0.23  ; OM2 V1
$OMEGA BLOCK(2)
 0.185  ; OM4 CL5
 0.178 0.258 ; OM5 CL4

$SIGMA
 0.371 ; SIG1 prop M
 0.206 ; SIG2 prop M3G
 0.0967 ; SIG3 prop M6G

$EST NOABORT PRINT=5 MAXEVAL=9999 METHOD=1 POSTHOC
$COV COMP
$TABLE ID TIME AMT RATE DV MDV CMT BWS PNA SGA CL2N CL2P CL3N CL3P CL4 CL5 V1 V2 V3 V4 Q ET1 ET2 ET3 ET4
NOAPPEND NOPRINT ONEHEADER FILE=sdtab5
