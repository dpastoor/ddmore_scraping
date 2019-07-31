$PROBLEM KPD model of CTC count and PSA

$INPUT ID TIME AMT DV CMT MDV

$DATA Simulated_KPD_CTC.count_PSA.csv IGNORE=@

$SUBROUTINE ADVAN13 TOL=9

$MODEL NCOMP=8
COMP=(A1) ;PK chemo
COMP=(A2) ;PK hormo
COMP=(A3) ;TS
COMP=(A4) ;CTC
COMP=(A5) ;Delayed PK chemo
COMP=(A6) ;Delayed PK hormo
COMP=(A7) ;Delayed TS
COMP=(A8) ;PSA

$PK CALLFL=-2   ; Call the PK subroutine with every event record, with additional and lagged doses
MU_1=LOG(THETA(1))
TS0=EXP(MU_1+ETA(1))      ; TS0

MU_2=LOG(THETA(2))
K1=EXP(MU_2+ETA(2))      ; Param PK Chemo

MU_3=LOG(THETA(3))
K2=EXP(MU_3+ETA(3))      ; Param PK Hormo

MU_4=LOG(THETA(4))
Q501=EXP(MU_4+ETA(4))      ; Param PK Chemo

MU_5=LOG(THETA(5))
Q502=EXP(MU_5+ETA(5))      ; Param PK Hormo

MU_6=LOG(THETA(6))
KOUTTS=EXP(MU_6+ETA(6))      ; Param TS

MU_7=THETA(7)
TH=EXP(MU_7+ETA(7))
KINTS=(TS0*KOUTTS)/((TH)/(1+(TH)))    ; Param TS

MU_8=LOG(THETA(8))
KINPSA=EXP(MU_8+ETA(8))

MU_9=LOG(THETA(9))
KOUTPSA=EXP(MU_9+ETA(9))

MU_10=LOG(THETA(10))
PSA0=EXP(MU_10+ETA(10))

MU_11=THETA(11)
K0=MU_11+ETA(11)    ; Zero-order production rate

MU_12=THETA(12)
ALAG5=MU_12+ETA(12)      ; TR: Delay duration
ALAG6=ALAG5
ALAG7=ALAG5

F4=K0*ALAG5              ; initial condition : R0=K0*TR . Bioaivalability for central compartment

MU_13=LOG(THETA(13))
OVDP=EXP(MU_13+ETA(13))

MU_14=LOG(THETA(14))
W1=EXP(MU_14+ETA(14))


A_0(1)=0
A_0(2)=0
A_0(3)=TS0
A_0(5)=0
A_0(6)=0
A_0(7)=TS0
A_0(8)=PSA0

$DES
DADT(1)=-K1*A(1)          ;time course of Chemo amount
DADT(2)=-K2*A(2)          ;time course of Hormo amount

DADT(3)=KINTS*(1-(A(1)/(Q501+A(1))))*(1-(A(2)/(Q502+A(2))))-KOUTTS*A(3)       ; Dynamic tumor size, latent variable

DADT(5)=-K1*A(5)           ; delayed time course of Chemo amount
DADT(6)=-K2*A(6)           ; delayed time course of Hormo amount

DADT(7)=KINTS*(1-(A(5)/(Q501+A(5))))*(1-(A(6)/(Q502+A(6))))-KOUTTS*A(7)        ; Delayed Dynamic tumor size, latent variable
A7=TS0
IF(T.GT.ALAG5) A7=A(7)

DADT(4)=K0*A(3)-K0*A7    ; time course of nber of cells
DADT(8)=KINPSA*A(3)-KOUTPSA*A(8)

$ERROR

NCTC=A(4)*0.0015       ; CTC in aliquots (alpha: scale factor)

PSA=A(8)
IF (PSA.LT.0.00001) PSA = 0.00001


CT=DV
IF(CT.LT.0) CT=0.00001


LFAC=GAMLN(CT+1.)

LGAM1=GAMLN(CT+1/OVDP)
LGAM2=GAMLN(1/OVDP)

LTRM1=(LOG(1/(1+OVDP*NCTC)))*(1/OVDP)
LTRM2=(LOG(NCTC/(NCTC+1/OVDP)))*(CT)

;Logarithm of the Negative Binomial distribution
LNB = LGAM1-LFAC-LGAM2+LTRM1+LTRM2  ;Ln(negative binomial)


IF (CMT.EQ.4) THEN
F_FLAG=2
;-2 Log Likelihood:
Y=-2*LNB
ENDIF


IF (CMT.EQ.8) THEN
F_FLAG=0
IPRED=LOG(PSA)
Y=IPRED+W1*ERR(1)
IRES=DV-IPRED
IWRES=IRES/W1
ENDIF



IF (ABS(ETA(1)).GT.50) EXIT 1 1
IF (ABS(ETA(2)).GT.50) EXIT 1 2
IF (ABS(ETA(3)).GT.50) EXIT 1 3
IF (ABS(ETA(4)).GT.50) EXIT 1 4
IF (ABS(ETA(5)).GT.50) EXIT 1 5
IF (ABS(ETA(6)).GT.50) EXIT 1 6
IF (ABS(ETA(7)).GT.50) EXIT 1 7
IF (ABS(ETA(8)).GT.50) EXIT 1 8
IF (ABS(ETA(9)).GT.50) EXIT 1 9
IF (ABS(ETA(10)).GT.50) EXIT 1 10
IF (ABS(ETA(11)).GT.100) EXIT 1 11
IF (ABS(ETA(12)).GT.20) EXIT 1 12
IF (ABS(ETA(13)).GT.50) EXIT 1 13


$THETA
1 FIX
2.48E-01  
4.49E-01  
2.61E-04  
3.97E-03  
5.13E-03  
6.33E+00  
1.40E+00  
8.13E-03  
1.53E+02  
3.08E+02  
5.77E+01
4.89E+00  
3.00E-01 
0.3

$OMEGA
0.0000001 FIX

$OMEGA BLOCK(9)
7.23E-01
-6.98E-01  1.83E+00
-6.38E-01 -6.34E-01  4.76E+00
-7.99E-01 -2.28E-01  2.93E+00  2.84E+00
3.51E-01  2.25E+00 -2.65E+00 -4.24E+00  2.06E+01
2.37E-01  4.00E-02 -7.51E-01 -9.54E-01  3.43E+00  7.86E-01
1.31E-01  9.74E-01 -1.81E+00 -1.45E+00 -6.04E-01 -9.75E-02  2.60E+00
3.12E-01  4.51E-02 -2.29E+00 -1.16E+00  1.84E+00  5.78E-01  9.67E-02  1.53E+00
1.15E-01  6.93E-01 -1.31E+00 -1.25E+00 -3.02E-01  1.53E-01  2.30E+00  3.48E-02  2.40E+00

$OMEGA BLOCK(2)
1.43E+03
2.94E+02  6.19E+01
$OMEGA
2.25

$OMEGA
0.0000001 FIX

$SIGMA
1 FIX ;ERR1
                                                                              
$ESTIMATION METHOD=SAEM LAPLACE INTER NUMERICAL SLOW NOHABORT NBURN=0 NITER=0 PRINT=1 NSIG=3 SIGL=9 GRD=DDDDDDDDDDDDDS
