$PROBLEM GENTAMICIN model
;data from Nielsen2009, Thomson1988, Estonia_unpub
$INPUT ID GA GIRL TIME RATE EVID AMT WT CREAT DV PNA PMA TCREA OCC
;GA in weeks
;AMT in mg
;RATE in mg/h
;CREAT in umol/L
;PNA in days
;PMA in weeks
;WT in g
;TCREA = typical (for PMA) SCr: TCREA=PMA*(-2.8488)+166.48 [Cuzzolin 2006 and Rudd 1983]
;OCC=a dose with a subsequent level reported
$DATA simdataDDM.csv IGNORE=@

$SUBROUTINE ADVAN6 TOL=6

$MODEL
      COMP=(CENTRAL)
      COMP=(PERIPH1)
      COMP=(PERIPH2)
      COMP=(COVCMT1)           ; PNA time-var. covariate compartment
      COMP=(COVCMT2)           ; CREATININE t-v. cov. compartment

$PK
; Three-comp model
IF(NEWIND.NE.2)OTIM1=0
IF(NEWIND.NE.2)OCOV1=0
IF(NEWIND.NE.2)OTIM2=0
IF(NEWIND.NE.2)OCOV2=0
;
STUDY=0
IF(ID.LT.2000) STUDY=1                     ;Glasgow, Thomson1988
IF(ID.GE.2000.AND.ID.LT.3000) STUDY=2      ;Uppsala, Nielsen2009
IF(ID.GE.3000) STUDY=3                     ;Estonia, unpublished
;
WTKG = WT/1000
;
T50  = THETA(7)
HILL = THETA(8)
MF   = PMA**HILL/(PMA**HILL+T50**HILL)
;
CREAT2 = CREAT
IF(CREAT.LT.0) CREAT2 = TCREA       ; when SCr is NA=-99, it is the typical SCr
;OF = (CREAT2/TCREA)**(THETA(9))
;
P50  = THETA(10)
;PNAF = PNA/(P50+PNA)
;
CRPWR = THETA(9)
;IOV code
BOVC = 0
IF(OCC.EQ.1)  BOVC = ETA(7)
IF(OCC.EQ.2)  BOVC = ETA(8)
IF(OCC.EQ.3)  BOVC = ETA(9)
IF(OCC.EQ.4)  BOVC = ETA(10)
IF(OCC.EQ.5)  BOVC = ETA(11)
IF(OCC.EQ.6)  BOVC = ETA(12)
IF(OCC.EQ.7)  BOVC = ETA(13)
IF(OCC.EQ.8)  BOVC = ETA(14)
IF(OCC.EQ.9)  BOVC = ETA(15)
IF(OCC.EQ.10) BOVC = ETA(16)
IF(OCC.EQ.11) BOVC = ETA(17)
IF(OCC.EQ.12) BOVC = ETA(18)
IF(OCC.EQ.13) BOVC = ETA(19)
IF(OCC.EQ.14) BOVC = ETA(20)
IF(OCC.EQ.15) BOVC = ETA(21)
IF(OCC.EQ.16) BOVC = ETA(22)
IF(OCC.EQ.17) BOVC = ETA(23)
IF(OCC.EQ.18) BOVC = ETA(24)
IF(OCC.EQ.19) BOVC = ETA(25)
IF(OCC.EQ.20) BOVC = ETA(26)
IF(OCC.EQ.21) BOVC = ETA(27)
IF(OCC.EQ.22) BOVC = ETA(28)
;
TVCL = THETA(1)*MF*(WTKG/70)**(0.632)   ; typical value of CL
TVV1 = THETA(2)*(WTKG/70)                       ; typical value of V1
TVQ  = THETA(3)*(WTKG/70)**(0.75)               ; ty. value of intercompartmental CL
TVV2 = THETA(4)*(WTKG/70)                       ; ty. value of V2
TVQ2 = THETA(5)*(WTKG/70)**(0.75)               ; ty value of CL3
TVV3 = THETA(6)*(WTKG/70)                       ; ty value of V3
;
CL   = TVCL*EXP(ETA(1)+BOVC)  ; individual value of CL
V1   = TVV1*EXP(ETA(2))
Q    = TVQ*EXP(ETA(3))
V2   = TVV2*EXP(ETA(4))
Q2   = TVQ2*EXP(ETA(5))
V3   = TVV3*EXP(ETA(6))
;
K    = CL/V1
K12  = Q/V1
K21  = Q/V2
K13  = Q2/V1
K31  = Q2/V3
;
IF(EVID.EQ.1) TM=TIME
IF(EVID.EQ.1) TAD=0
IF(EVID.NE.1) TAD=TIME-TM
;
SL1 = 0
IF(TIME.GT.OTIM1) SL1 = (PNA-OCOV1)/(TIME-OTIM1)
A_0(4) = PNA
;
SL2 = 0
IF(TIME.GT.OTIM2) SL2 = (CREAT2-OCOV2)/(TIME-OTIM2)
A_0(5) = CREAT2
;
$DES
DADT(4)= SL1
TCOV1 = A(4)

DADT(5)= SL2
TCOV2 = A(5)

;PNAF = PNA/(P50+PNA)
PNAF = TCOV1/(P50+TCOV1)
;OF = (CREAT2/TCREA)**(THETA(9))
OF = (TCOV2/TCREA)**CRPWR
DADT(1) = A(3)*K31+A(2)*K21-A(1)*(K*PNAF*OF+K12+K13)
DADT(2) = A(1)*K12-A(2)*K21
DADT(3) = A(1)*K13-A(3)*K31

$ERROR
 IPRED  = A(1)/V1
 Y      = IPRED*(1+EPS(1)) + EPS(2)

OCOV1 = PNA
OTIM1 = TIME

OCOV2 = CREAT2
OTIM2 = TIME

$THETA  (0,6.20684) ; 1. TVCL (lower bound,initial estimate)
$THETA  (0,26.5004) ; 2. TVV1  (lower bound,initial estimate)
$THETA  (0,2.15099) ; 3. TVQ
$THETA  (0,21.151) ; 4. TVV2
$THETA  (0,0.270697) ; 5. TVQ2
$THETA  (0,147.893) ; 6. TVV3
$THETA  55.4 FIX ; 7. T50
$THETA  3.33 FIX ; 8. Hill
$THETA  -0.129934 ; 9. power exponent on creatinine
$THETA  (0,1.70302) ; 10. PNA50
$OMEGA  BLOCK(2)
 0.175278  ; variance for ETA(1), initial estimate
 0.115896 0.112362  ; COvariance ETA(1)-ETA(2), var for ETA(2), initial estimate
$OMEGA  0  FIX
$OMEGA  0.131759
$OMEGA  0  FIX
$OMEGA  0.177214
$OMEGA  BLOCK(1)
 0.0140684  ;  7. IOV_CL
$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$OMEGA  BLOCK(1) SAME

$SIGMA  0.036033  ; variance PROP res error, initial estimate
$SIGMA  0.0164023
$ESTIMATION METHOD=1 INTER MAXEVAL=0 PRINT=1  ; calculation method
;$COVARIANCE                                      ; standard error of estimate is calculated

$TABLE ID TIME IPRED DV CWRES 
       CL V1 Q V2 Q2 V3 ETA(1) ETA(2) ETA(3) ETA(4) ETA(5) ETA(6)
       GA GIRL RATE EVID AMT WT CREAT PNA PMA TCREA OCC TAD STUDY 
       MF OF CREAT2 PNAF BOVC SL2 TCOV2 PNA SL1 TCOV1
       TVCL TVV1 TVQ TVV2 TVQ2 TVV3                             NOPRINT ONEHEADER FILE=sdtab35b_ddm2
;$TABLE ID WT GA PNA PMA CREAT TCREA                             NOPRINT NOAPPEND ONEHEADER FILE=cotab32
;$TABLE ID GIRL STUDY                                            NOPRINT NOAPPEND ONEHEADER FILE=catab32
;$TABLE ID CL V ETA(1) ETA(2)       NOPRINT NOAPPEND ONEHEADER FILE=patab27
;$SCAT DV VS PRED UNIT

