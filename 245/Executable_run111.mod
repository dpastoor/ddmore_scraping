;Rerun to estimate parameters at CLCR 70 and WT 75
;;1. Aim with run
;;   Based on run93
;;   Adding BLOCK(2)
;;2. Based on run31 from NXY-0004
;;3. Structural model
;;   Two-compartment model
;;4. Covariate model
;;   CLCR on CL (nonlinear)
;;   WT on V1
;;5. Interindividual model
;;   Exponential IIV for CL, V1,
;;   and correlation between CL and V1
;;6. Residual variability
;;   Proportional (but add w logtr data)
;;7. Estimation
;;   FOCE
;;8. Other
;;   Simulated_comb2.dta
$PROB NXY-059 Analysis, SA-NXY-0004 and 0003
$INPUT ID OID=DROP CENT=DROP PNO=DROP TARG DAT2=DROP 
       TIME AMT DUR=DROP RATE ODV DV FU FLAG=DROP
       SEX AGE WT HT BMI RACE SCR CLCR FLA2 STUD
$DATA /home/pm_common/Model_database/DDMoRe_Prep/Jonsson-S-2005_NXY059/Data/Simulated_comb2.dta IGNORE=@ 
$SUBROUTINES ADVAN3 TRANS4
$PK
 CREA = CLCR
 IF(CLCR.EQ.-99) CREA = 61.48
 IF(CLCR.EQ.-99.AND.ID.EQ.10) CREA = 34.55
 IF(CLCR.EQ.-99.AND.ID.EQ.21) CREA = 124.22
 IF(CREA.LE.40)  CLCLCR = 0
 IF(CREA.GT.40)  CLCLCR = THETA(6)*(CREA-40)

 IF(WT.EQ.-99) THEN
    V1WT = 0
 ELSE
    V1WT = THETA(7)*(WT-76.00)
 ENDIF

CLCOV=1+CLCLCR
V1COV=1+V1WT

TVCL   = THETA(5)*CLCOV
TVV1   = THETA(2)*V1COV
TVQ    = THETA(3)
TVV2   = THETA(4)

CL     = TVCL*EXP(ETA(1))
V1     = TVV1*EXP(ETA(2))
Q      = TVQ
V2     = TVV2

S1     = V1

$ERROR
NPRE   = F+0.000001
IPRED  = LOG(F+0.000001)
IRES   = DV-IPRED
W      = THETA(1)
IWRES  = IRES/W
Y      = IPRED+EPS(1)*W


$THETA (0,0.1650)                   ;1 prop error
$THETA (0,7.8500)                   ;2 V1
$THETA (0,13.100)                   ;3 Q
$THETA (0,7.2000)                   ;4 V2
$THETA (0,2.8800)                   ;5 TVCL
$THETA (-0.00972,0.0192,0.0505)     ;6 CLCR on CLR
$THETA (-0.020,0.0184,0.027)        ;7 WT on V1       

$OMEGA BLOCK(2) 0.0536 0.02 0.160   ;IIV CL, V1

$SIGMA 1 FIX
$EST PRINT=2 METHOD=1 SIG=3 MSFO=msfb111 NOABORT MAX=0
;;COV
;;$TAB ID TIME TARG STUD ODV NPRE IPRED IWRES  ONEHEADER NOPRINT FILE=sdtab111
;;$TAB ID V1 ETA(1) CL ETA(2) TVCL TVV1        ONEHEADER NOPRINT FILE=patab111
;;$TAB ID FU AGE WT HT BMI SCR CLCR            ONEHEADER NOPRINT FILE=cotab111
;;$TAB ID SEX RACE                             ONEHEADER NOPRINT FILE=catab111
