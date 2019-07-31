;; 1. Based on: run06
;; 2. Description: Paracetamol (APAP) in term and preterm newborns
;; x1. Author: user
$PROBLEM Paracetamol (APAP) in term and preterm newborns
$INPUT 
ID 
TIME ; in min
AMT  ; in mg
RATE ; in mg/min
DV   ; in blood mg/L in urine mg
CMT  ; see under $MODEL and in publication
MDV 
EVID 
BWS  ; bodyweight at time of study in kg

$DATA Simulated_ParacetamolPKnewborns.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=5

$MODEL
NCOMPARTMENTS=6
COMP=(CENTRAL,DEFDOSE) ;PLASMA APAP
COMP=                  ;PLASMA APAP-G
COMP=                  ;PLASMA APAP-S
COMP=                  ;URINE APAP-G
COMP=                  ;URINE APAP
COMP=                  ;URINE APAP-S

$PK
ET1=EXP(ETA(1))
ET2=EXP(ETA(2))
ET3=EXP(ETA(3))
ET4=EXP(ETA(4))

TH1=THETA(1)
TH2=THETA(2)/1000
TH3=THETA(3)/1000
TH4=THETA(4)/1000
TH5=THETA(5)
TH6=THETA(6)

V1=TH1*BWS*ET1         ;volume plasma APAP
V2=V1*0.18             ;volume metabolites
V3=V2

CLG=TH2*BWS*ET2        ;formation clearance APAP-G
CLS=TH3*BWS**TH6*ET3   ;formation clearance APAP-S
CLA=TH4*BWS*ET4        ;clearance unchanged APAP 

K12=CLG/V1
K13=CLS/V1
K15=CLA/V1
K24=TH5*K15
K36=K24

CLG2=K24*V2            ;elimination clearance APAP-G
CLS2=K36*V3            ;elimination clearance APAP-S

S1=V1
S2=V2
S3=V3
S4=1
S5=1
S6=1

$DES
DADT(1)=-K12*A(1)-K13*A(1)-K15*A(1)
DADT(2)=K12*A(1)-K24*A(2)
DADT(3)=K13*A(1)-K36*A(3)
DADT(4)=K24*A(2)
DADT(5)=K15*A(1)
DADT(6)=K36*A(3)

$ERROR
COM1=0
IF(CMT.EQ.1) COM1=1
COM4=0
IF (CMT.EQ.4) COM4=1
COM5=0
IF (CMT.EQ.5) COM5=1
COM6=0
IF (CMT.EQ.6) COM6=1
IPRE=F

Y1=F*(1+ERR(1))+ERR(5) ;Y1=APAP PLASMA
Y4=F*(1+ERR(2))        ;Y4=APAP-G URINE
Y5=F*(1+ERR(3))        ;Y5=APAP URINE
Y6=F*(1+ERR(4))        ;Y6=APAP-S URINE

Y=COM1*Y1+COM4*Y4+COM5*Y5+COM6*Y6

$THETA
(0, 1.06) ; TH1 V1 APAP
(0, 0.266) ; TH2 CL APAP-G
(0, 1.46) ; TH3 CL APAP-S
(0, 0.285) ; TH4 CL APAP
(0, 11.3) ; TH5 CL 24 EN 36
(0, 1.4) ; TH6 EXP on CLS

$OMEGA
 0.0925 ; V
 0.599 ; CLAPAP-G
 0.312 ; CLAPAP-S
 0.0879 ; CLAPAP

$SIGMA
 0.0198 ; 1 APAP plasma proportional
 0.223 ; 2 APAP-G URINE proportional
 0.188 ; 3 APAP URINE proportional
 0.332 ; 4 APAP-S URINE proportional
 0.354 ; 5 APAP PLASMA additive

$ESTIMATION PRINT=1 MAXEVAL=9999 METHOD=CONDITIONAL INTERACTION NOABORT
$COV PRINT=E

$TABLE ID TIME AMT RATE DV CMT MDV EVID BWS IPRE NOAPPEND ONEHEADER NOPRINT FILE=sdtab07
