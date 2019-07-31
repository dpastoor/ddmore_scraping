$PROBLEM + COVAR KM-VMAX
; 1. Based on: run530
; 2. Description:
;    DOSE ON F AS EMAX MODEL
; 3. Label:
;    HIGHRIF1 PK MODEL
; 4. Structural model:
;    1-COMP PK MODEL WITH TRANSIT ABSORPTION MODEL WITH AUTO-INDUCTION AND MICHAELIS-MENTEN CL WITH M3 METHOD LOG-TRANSFORMED DATA
; 5. Covariate model:
;    FFM ON CL AND V2
; 6. Inter-individual variability:
;    KM, V2, MTT, NN, KA, VMAX
; 7. Inter-occasion variability:
;    KM, MTT, BIO, KA.
; 8. Residual variability:
;    ADDITIVE
; 9. Estimation:
;    LAPLACIAN WITH INTERACTION

; Based on final model by Smythe run 106

$ABBREVIATED DERIV2=NOCOMMON

$INPUT ID TIME TADO DGRP DV BQL AMT EVID OCC PLOT AGE SEX RACE WT HT BMI HIV FFM NDV DOSE
; TIME=hours, TADO=time after last dose, DGRP=dose in mg/kg,
; BQL (0=observation is not BLOQ, 1=observation is BLOQ, 2=observation missing, 3=dummy or dosing time point)
; OCC (1=day 7, 2=day 14), PLOT=flag variable for creating VPCs, SEX (1=male, 0=female)
; FFM=fat-free mass in kg, DOSE=dose in mg
$DATA       Simulated_Rif_PK_data.csv IGNORE=@
$SUBROUTINE ADVAN13 TRANS1 TOL=10
$MODEL      NCOMP = 3 COMP = (DEPOT,DEFDOSE) COMP = (CENTRAL,DEFOBS) COMP = (ENZ)

$PK
"FIRST
"  COMMON/PRCOMG/  IDUM1,IDUM2,IMAX,IDUM4,IDUM5
"  INTEGER IDUM1,IDUM2,IMAX,IDUM4,IDUM5
"  IMAX=1000000000

IF(AMT.GT.0)PD     = AMT          ; PD   = oral DOSE
IF(AMT.GT.0)TDOS   = TIME         ; TDOS = time of DOSE
TAD                = TIME - TDOS  ; TAD  = time after DOSE

;;-------------------Fat Free Mass-------------------------------------;;

NFMCL  = FFM
ALLMCL = (NFMCL/70)**0.75

NFMV  = FFM
ALLMV = (NFMV/70)

TVVMAX  = THETA(1)               ; CL
TVKM    = THETA(2)
TVV2    = THETA(3)*ALLMV               ; V2
TVKA    = THETA(4)                     ; KA
TVEMAX  = THETA(5)                     ; EMAX
TVEC50  = THETA(6)                     ; EC50
TVKENZ  = THETA(7)                     ; KENZ
TVMTT   = THETA(8)                     ; MTT
TVNN    = THETA(9)                     ; NN
TVFEMAX = THETA(10)                    ; FEMAX
TVFED50 = THETA(11)                    ; FEC50

IF (OCC.EQ.1) THEN
  IOVBIO = ETA(7)
  ELSE
  IOVBIO = ETA(8)
ENDIF

IF (OCC.EQ.1) THEN
  IOVMTT = ETA(9)
  ELSE
  IOVMTT = ETA(10)
ENDIF

IF (OCC.EQ.1) THEN
  IOVKM = ETA(11)
  ELSE
  IOVKM = ETA(12)
ENDIF

IF (OCC.EQ.1) THEN
  IOVKA = ETA(13)
  ELSE
  IOVKA = ETA(14)
ENDIF

VMAX    = TVVMAX*EXP(ETA(2))
KM      = TVKM*EXP(ETA(1)+IOVKM)
V2      = TVV2*EXP(ETA(3))
KA      = TVKA*EXP(ETA(6)+IOVKA)
EC50    = TVEC50
EMAX    = TVEMAX
KENZ    = TVKENZ
FEMAX   = TVFEMAX
FED50   = TVFED50
TVBIO   = 1*(1+FEMAX*(DOSE-450)/(FED50+(DOSE-450)))
BIO     = TVBIO*EXP(IOVBIO)
K       = CL/V2
MTT     = TVMTT*EXP(ETA(4)+IOVMTT)
NN      = TVNN*EXP(ETA(5))
S2      = V2

F1      = 0 ; Transit absorption compartment
A_0(2)  = 0.0001 ; Central comp
A_0(3)  = 1 ; Induction compartment

KTR     = (NN + 1) / MTT

L       = 0.9189385 + (NN + 0.5)*LOG(NN) - NN + LOG(1 + 1/(12*NN)) ; logarithm of the approximation to the gamma function
LBPD  = LOG(BIO*PD)
LKTR  = LOG(KTR)
CUMUL = LBPD + LKTR - L

$DES
CP      = A(2)/V2

TEMPO   = T - TDOS
  IF(TEMPO.GT.0)THEN
	KTT   = KTR*TEMPO
  DADT(1) = EXP(CUMUL + NN*LOG(KTT) - KTT) - KA*A(1)
  ELSE
        KTT   = 0
  DADT(1) = 0
  ENDIF

DADT(2) = KA*A(1) - (((VMAX/(KM+CP))*ALLMCL)/V2)*A(2)*A(3)
EFF     = (EMAX*(CP)) / (EC50 + CP)
DADT(3) = KENZ*(1 + EFF) - KENZ*A(3)

$ERROR
IPRED   = LOG(A(2)/S2+0.00001)
ADD     = SQRT(SIGMA(1,1))   ; ADD error
SD      = SQRT((ADD)**2)

;Sim_start
LLOQ=LOG(0.13)

DUM=(LLOQ-IPRED)/SD
CUMD=PHI(DUM)

IF(DV.GE.LLOQ) THEN
F_FLAG = 0
IRES   = DV - IPRED
IWRES  = IRES / SD
Y      = IPRED + EPS(1)
ELSE
F_FLAG = 1
IRES   = 0
IWRES  = 0
MDVRES = 1
Y=CUMD
ENDIF

;IRES  = DV - IPRED
;IWRES = IRES / SD
;Y     = IPRED + EPS(1)
;Sim_end

AA1 = A(1)     ; absorption comp
AA2 = A(2)     ; central rif comp
AA3 = A(3)     ; auto-induction comp

IF(AMT.GT.0) THEN
	TDOS = TIME
	PD   = AMT
ENDIF

$THETA  (0,525) ; 1 VMAX
$THETA  (0,35.3) ; 2 KM
$THETA  (0,87.2) ; 3 V2
$THETA  (0,1.77) ; 4 KA
$THETA  (0,1.16,1.25) ; 5 EMAX
$THETA  (0,0.0699) ; 6 EC50
$THETA  (0.005,0.00603) ; 7 KENZ
$THETA  (0,0.513) ; 8 MTT
$THETA  (1,23.8) ; 9 NN
$THETA  (0,0.504) ; 10 FEMAX
$THETA  (0,67) ; 11 FED50
$OMEGA  BLOCK(2)
 0.128  ; 1 IIV in KM
 0.0418 0.0901  ; 2 IIV IN VMAX
$OMEGA  0.00618  ; 3 IIV in V2
$OMEGA  0.146  ; 4 IIV in MTT
$OMEGA  0.607  ; 5 IIV in NN
$OMEGA  0.114  ; 6 IIV in KA
$OMEGA  BLOCK(1)
 0.0248  ; 7 IOV in F
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1)
 0.318  ; 9 IOV in MTT
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1)
 0.0355  ; 11 IOV in KM
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1)
 0.0985  ; 13 IOV in KA
$OMEGA  BLOCK(1) SAME
$SIGMA  0.0555  ;  ADD ERROR

;Sim_start
$ESTIMATION METHOD=1 LAPLACIAN INTER NUMERICAL SLOW MAXEVAL=0 NSIG=3 SIGL=9 PRINT=3 MCETA=100
$COVARIANCE PRINT=E MATRIX=S SLOW
;$SIMULATION (1234) ONLYSIM
;Sim_end
$TABLE ID IPRED IWRES CWRES NPDE DV OCC TIME TADO DGRP PLOT CP
       NOPRINT ONEHEADER FILE=sdtab537
$TABLE ID
       NOPRINT ONEHEADER FILE=catab537
$TABLE ID AGE SEX RACE WT HT BMI HIV FFM
       NOPRINT ONEHEADER FILE=cotab537
$TABLE ID V2 MTT BIO NN KM KA VMAX ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7 ETA8 ETA9 ETA10 ETA11 ETA12 ETA13
       NOPRINT ONEHEADER FILE=patab537
$TABLE ID
       NOPRINT ONEHEADER FILE=mytab537