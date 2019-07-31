$SIZES      NO=800 LIM6=800     ; needed for big datasets
$PROB Tumor model

; Reference: Hansson E.K. 2011. Pharmacometric models for Biomarkers, Side Effects and Efficacy in Anticancer Drug Therapy. Acta Universitatis Upsaliensis.

$INPUT
ID                    ; Patient identification
CYCLE                 ; Cycle number
TIME                  ; Time in hours
WEEK                  ; Time in weeks
FLAG                  ; 1. Dose ; 4. Tumor size (SLD)
DV
DOS                   ; Sunitinib dose
PLA                   ; Placebo: 1. untreated, 0. treated
CL                    ; Posthoc total plasma clearance
EVID                  ; 0. observation; 2. other event
BAS3 MRT3 EC53        ; posthoc parameters for VEGFR3 timecourse
SBAS SMRT SEC5 SLO    ; posthoc parameters for SKIT timecourse



$DATA   data_ddmore_TGI_GIST_simulated.csv
        IGN=#

$SUBROUTINE ADVAN6 TOL=5
$MODEL  NCOMP=4

$PK

; Verbatim code : changes the iteration maximum (IMAX) (default value 100000)
"FIRST
" COMMON /PRCOMG/ IDUM1,IDUM2,IMAX,IDUM4,IDUM5
" INTEGER IDUM1,IDUM2,IMAX,IDUM4,IDUM5
" IMAX=100000000

;-----DRUG EXPOSURE------------------
  AUC=DOS/CL

;-----sKIT---------------------------
  BM0S    =    SBAS                       ; model predicted parameters for SKIT timecourse
  MRTS    =    SMRT
  IMAX1    =    1
  IC50S   =    SEC5
  DPSLOS  =    SLO
  KOUTS   =    1/MRTS

;-----VEGFR3-------------------------
  BM03    =    BAS3                       ; model predicted parameters for SKIT timecourse
  IC503   =    EC53
  KOUT3   =    1/MRT3


 IF(TIME.EQ.0.AND.FLAG.EQ.4)THEN
   OBASE  =    DV                         ; observed tumor size at baseline (time 0)
 ENDIF

 W1       =    THETA(4)*OBASE
 IBASE    =    OBASE+ETA(5)*W1            ; observed tumor size at baseline acknowledging residual error

  TVKG    =    THETA(1)/24/7
  KG      =    TVKG  *EXP(ETA(1))         ; tumor growth rate constant
  TVKSKI  =    THETA(2)/24/7
  KSKIT   =    TVKSKI*EXP(ETA(2))         ; tumor size reduction rate constant related to SKIT response
  TVLAM   =    THETA(3)/24/7
  LAMBDA  =    TVLAM *EXP(ETA(3))         ; resistance appearance rate constant
  TVDRU   =    THETA(5)/24/7
  KDRUG   =    TVDRU *EXP(ETA(4))         ; exposure driven drug effect
  TVKV3   =    THETA(6)/24/7
  KVEG3   =    TVKV3                      ; tumor size reduction rate constant related to sVEGFR3 response

; NB: THETAs are divided by 24 and 7 to scale the rate constants from 1/week to 1/hour.

;-----Compartment initialization-----
A_0(1) = BM0S                             ; SKIT
A_0(2) = BM0S                             ; SKIT
A_0(3) = BM03                             ; VEGFR3
A_0(4) = IBASE                            ; TUMOR


$DES
;-----SKIT---------------------------
  EFFS    =    IMAX1*AUC/(IC50S+AUC)       ; inhibitory Emax drug effect on SKIT
  DPS     =    BM0S*(1+DPSLOS*T)          ; time-dependent linear disease progression model describing the change in SKIT in untreated patients
  KINS    =    DPS*KOUTS

;-----VEGFR3-------------------------
  EFF3    =    IMAX1*AUC/(IC503+AUC)       ; inhibitory Emax drug effect on VEGFR3
  KIN3    =    BM03*KOUT3


DADT(1)   =    KINS*(1-EFFS)-KOUTS*A(1)   ; SKIT timecourse with drug effect
DADT(2)   =    KINS-KOUTS*A(2)            ; SKIT timecourse in untreated patients
DADT(3)   =    KIN3*(1-EFF3)-KOUT3*A(3)   ; VEGFR3 timecourse


;-----TUMOR--------------------------
  SKIT    =    ((A(1)-A(2))/A(2))*KSKIT   ; effect of SKIT (SKIT timecourse relative to baseline)
  VEG3    =    ((A(3)-BM03)/BM03)*KVEG3 ; effect of VEGFR3 (VEGFR3 timecourse relative to baseline)

  AUC1= AUC*KDRUG                         ; exposure driven drug effect


DADT(4)   =    KG*A(4)-(AUC1+(-SKIT)+(-VEG3))*EXP(-(LAMBDA*T))*A(4)   ; longitudinal model describing tumor growth

$ERROR
  AA1 = A(1)                              ; SKIT
  AA2 = A(2)                              ; SKIT
  AA3 = A(3)                              ; VEGFR3
  TUM = A(4)                              ; TUMOR


STRT=0                                    ; stratification variable
IF(PLA.EQ.1)STRT=1



  IF(FLAG.EQ.4) THEN
  IPRED   =    A(4)
  W       =    IPRED*THETA(4)
  Y       =    IPRED+W*EPS(1)             ; proportional residual error
  IRES    =    DV -  IPRED
  IWRES   =    IRES/W
  ENDIF

;-----INITIAL ESTIMATES------------------------

$THETA (0,0.0118)                         ; 1 KG
$THETA (0,0.00282                         ; 2 KD SKIT
$THETA (0,0.0217)                         ; 3 LAMBDA
$THETA (0,0.125)                          ; 4 RES
$THETA (0,0.00503)                        ; 5 KDRUG
$THETA (0,0.0371)                         ; 6 KD VEGFR3


$OMEGA 0.29                               ; 1 KG
$OMEGA 5.91                               ; 2 KD SKIT
$OMEGA 0 FIX                              ; 3 LAMBDA
$OMEGA 1.42                               ; 4 KDRUG
$OMEGA 1 FIX                              ; 5 IBASE

$SIGMA 1 FIX

$EST PRINT=1 MAXEVAL=0 METHOD=1 INTER POSTHOC MSFO=MSF2        ; FOCE with interaction
;$COV
$TABLE ID TIME IPRED IBASE IWRES IRES AUC LAMBDA KG KSKIT KVEG3 AA1 AA2 AA3 TUM  ETA1 ETA2 NOPRINT ONEHEADER FILE=sdtab2