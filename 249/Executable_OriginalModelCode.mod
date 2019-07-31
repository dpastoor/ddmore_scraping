$PROBLEM Midazolam PK in critically ill Children

$INPUT 
ID 
TIME ; in hours
CMT 
AMT  ; in microgram
RATE ; in microgram / hour
DV   ; in microgram / L
MDV 
ORGF ; number of organs that are failig 
WT   ; bodyweight in kg
CRP  ; CRP concentration in mg/L

$DATA Data20150701.csv IGNORE=@

$SUBROUTINES ADVAN6 TOL=6 

$MODEL
COMP=(CENTRAL)
COMP=(PERIPH)

$PK
OCC=1
IF(TIME.GE.24)OCC=2
IF(TIME.GE.48)OCC=3
IF(TIME.GE.72)OCC=4
IF(TIME.GE.96)OCC=5
IF(TIME.GE.120)OCC=6

FLAG1=0 
FLAG2=0 
FLAG3=0 
FLAG4=0 
FLAG5=0 
FLAG6=0 

IF(OCC.EQ.1) FLAG1 = 1
IF(OCC.EQ.2) FLAG2 = 1
IF(OCC.EQ.3) FLAG3 = 1
IF(OCC.EQ.4) FLAG4 = 1
IF(OCC.EQ.5) FLAG5 = 1
IF(OCC.EQ.6) FLAG6 = 1

IOV=FLAG1*ETA(3)+FLAG2*ETA(4)+FLAG3*ETA(5)+FLAG4*ETA(6)+FLAG5*ETA(7)+FLAG6*ETA(8)

IIVCL = ETA(1)

ETCL = IIVCL + IOV

IF (ORGF.EQ.0) TVCL = THETA(1) * (WT/5)**(THETA(5)) * (CRP/32)**THETA(11)
IF (ORGF.EQ.1) TVCL = THETA(7) * (WT/5)**(THETA(5)) * (CRP/32)**THETA(11)
IF (ORGF.EQ.2) TVCL = THETA(8) * (WT/5)**(THETA(5)) * (CRP/32)**THETA(11)
IF (ORGF.EQ.3) TVCL = THETA(9) * (WT/5)**(THETA(5)) * (CRP/32)**THETA(11)
IF (ORGF.GT.3.5) TVCL = THETA(10) * (WT/5)**(THETA(5)) * (CRP/32)**THETA(11)
CL = TVCL * EXP(ETCL)
TVV1 = THETA(2) * (WT/5)**(THETA(6))
V1 = TVV1 * EXP(ETA(2))
TVQ  = THETA(3) 
Q  = TVQ  
TVV2 = THETA(4) 
V2 = TVV2 
S1 = V1

K10 = CL/V1
K12 = Q/V1
K21 = Q/V2

$DES
DADT(1) = -K10*A(1) - K12*A(1) + K21*A(2)
DADT(2) =             K12*A(1) - K21*A(2)

$ERROR
Y=F*(1+ERR(1))+ERR(2)
IPRED=F

$THETA
1.6 FIX       ;CL orgf = 0
(0.001, 10 )  ;V1 (L)
(0.001,  1 )  ;Q (L/h)
(0.001, 10 )  ;V2 (L)
(0.001, 0.75) ;WT-CL pow
(0.001, 0.75) ;WT-V1 pow
(0.001,  1 )  ;CL orgf 1
(0.001,  1 )  ;CL orgf 2
(0.001,  1 )  ;CL orgf 3
(0.001,  1 )  ;CL orgf 4/5
(-5, -0.5)    ;CRP-CL pwr

$OMEGA
0.1 ;IIV CL
0.1 ;IIV V1
$OMEGA BLOCK(1) 0.05 ;IOV 
$OMEGA BLOCK(1) SAME ;IOV day2		
$OMEGA BLOCK(1) SAME ;IOV day3		
$OMEGA BLOCK(1) SAME ;IOV day4		
$OMEGA BLOCK(1) SAME ;IOV day5		
$OMEGA BLOCK(1) SAME ;IOV day6

$SIGMA
0.1 ;Proportional error PK
0.1 ;Additive error PK

$EST METHOD=1 INTER MAXEVAL=2000 NOABORT SIG=3 PRINT=5 POSTHOC
$COV PRINT=E

$TABLE ID TIME DV MDV PRED IPRED CWRES WT ONEHEADER NOPRINT FILE=sdtab
$TABLE ID CL V1 V2 Q ETA(1) ETA(2) ETCL IOV WT ORGF CRP ONEHEADER NOPRINT FIRSTONLY FILE=patab


