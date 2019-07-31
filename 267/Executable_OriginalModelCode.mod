$PROB Propacetamol (APAP) in young (pregnant) women
$INPUT 
ID 
TIME ; in hours
MDV 
RATE ; in mg/h
AMT  ; in mg
CMT  ; see definitions under $MODEL
EVID 
DV   ; mg/L (plasma) or mg (urine)
BW   ; in kg
TERM ; term birth yes(1) / no(0)
OCC  ; 1=pregnant, 2=2 weeks post delivery, 3=1 year post delivery, 4=volunteer no BC, 5=volunteer on BC
BC   ; on Birth control yes(1) / no(0)
UF   ; urine flow mL/h

; CL and Q = L/h
; V = L
; ka = /h

$DATA Simulated_APAP_YoungWomen.csv IGNORE=@ 
$SUBROUTINE ADVAN6 TOL=5

$MODEL
NCOMPARTMENTS=8
COMP=(CENTRAL,DEFDOSE) ;1 Plasma central APAP
COMP=(PERIAPAP)        ;2 Plasma perpheral APAP 1
COMP=(PAPAPG)          ;3 PLASMA APAP-Glucuronide
COMP=(PAPAPS)          ;4 PLASMA APAP-Sulphate
COMP=(UAPAPG)          ;5 URINE APAP-Glucuronide
COMP=(UAPAPS)          ;6 URINE APAP-Sulphate
COMP=(UAPAP)           ;7 URINE APAP
COMP=(PERI2)           ;8 Plasma perpheral APAP 2 

$PK
V2=THETA(5)       
V1=THETA(8)*EXP(ETA(1))                 
IF(OCC.EQ.1) V1=(THETA(1)*THETA(8))*EXP(ETA(1))
Q1=(THETA(6)*((BW/70)**1)) 
V3=V1*0.18       
V4=V3
NEO=0
IF(BC.EQ.1) NEO=1        
CLG=THETA(9)*EXP(ETA(2))
IF(OCC.EQ.1) CLG=(THETA(2)*THETA(9))*EXP(ETA(2))
IF(OCC.EQ.2) CLG=THETA(11)*THETA(9)*EXP(ETA(2))
ICLG=CLG
IF(BC.EQ.1) CLG=THETA(12)*ICLG
CLS=(TERM*THETA(3)+(1-TERM)*THETA(10))*EXP(ETA(3))
RCUF=(THETA(16)*(UF-100))
IF(UF.EQ.0) RCUF=0
CLA=(THETA(4)+RCUF)*EXP(ETA(4))
V8=THETA(13)
Q2=THETA(14)
IF(OCC.EQ.2) Q2=(THETA(15)*THETA(14))
CLT=CLS+CLG+CLA

K12=Q1/V1
K21=Q1/V2
K18=Q2/V1
K81=Q2/V8
K13=CLG/V1
K14=CLS/V1
K17=CLA/V1
K35=THETA(7)*K17
K46=K35
CLG2=K35*V3     
CLS2=K46*V4    

S1=V1
S3=V3
S4=V4
S5=1
S6=1
S7=1

$DES
DADT(1)=-K12*A(1)+K21*A(2)-K13*A(1)-K14*A(1)-K17*A(1)-K18*A(1)+K81*A(8)
DADT(2)=K12*A(1)-K21*A(2)
DADT(3)=K13*A(1)-K35*A(3)
DADT(4)=K14*A(1)-K46*A(4)
DADT(5)=K35*A(3)
DADT(6)=K46*A(4)
DADT(7)=K17*A(1)
DADT(8)=K18*A(1)-K81*A(8)

$ERROR
COM1=0
IF(CMT.EQ.1) COM1=1
COM5=0
IF(CMT.EQ.5) COM5=1
COM6=0
IF(CMT.EQ.6) COM6=1
COM7=0
IF(CMT.EQ.7) COM7=1
GREG=0
IF(OCC.EQ.5) GREG=1
IPRE=F

Y1=F*(1+ERR(1))
Y2=F*(1+ERR(5))+ERR(6)
Y5=F*(1+ERR(2))        
Y6=F*(1+ERR(3))       
Y7=F*(1+ERR(4))      
Y=(1-GREG)*(COM1*Y1+COM5*Y5+COM6*Y6+COM7*Y7)+ GREG*Y2

$THETA
( 0,   1.86   )   ;fractional V1 difference OCC=1 TH1
( 0,   2.03   )   ;fractional CL formation APAP-G for OCC=1 TH2
( 0,   5.61   )   ;CL formation APAP-S for TERM=0 TH3
( 0,   0.925  )   ;CL of unchanged APAP TH4
( 0,  19.7    )   ;V2 TH5
( 0,   1.29   )   ;Q1 TH6
( 0,   4.62   )   ;k35 as fraction of k17 TH7
( 0,  18.5    )   ;V1 for OCC > 1TH8
( 0,   7.33   )   ;CL formation APAP-G TH9
( 0,   3.86   )   ;CL formation APAP-S for TERM=1 TH10
( 0,   0.547  )   ;fractional CL formation APAP-G for OCC=2 TH11
( 0,   1.46   )   ;fractional CL formation APAP-G for BC=1 TH12
( 0,  23.9    )   ;V8 TH13
( 0,  61.1    )   ;Q2 if OCC!=2 TH14
( 0,   0.128  )   ;fractional Q2 for OCC=2 TH15
( 0,   0.00535)   ;slope impact UF TH16

$OMEGA 
0.0867    ;V1 ETA1
0.121     ;CL fomation APAP-G ETA2
0 FIXED   ;CL fomation APAP-S ETA3
0.122     ;CL unchanged APAP ETA4

$SIGMA
0.0695  ;APAP plasma PROP EPS1
0.292   ;APAP-G URINE PROP EPS2
0.147   ;APAP-S URINE PROP EPS3
0.152   ;APAP URINE PROP EPS4
0.0169  ;PROP error in OCC=5 EPS5
0.016   ;ADD error in OCC=5 EPS6

$EST NOABORT SIGDIG=3 PRINT=5 MAXEVAL=9999 METHOD=1 INTERACTION POSTHOC
$COV COMP
$TABLE ID TIME MDV RATE AMT CMT EVID DV BW TERM OCC BC UF V1 V2 V8 Q1 Q2 CLG CLS CLA CLT IPRE PRED FILE=sdtabSimulatedData NOPRINT NOAPPEND

