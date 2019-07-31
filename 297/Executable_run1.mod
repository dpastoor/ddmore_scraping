;; 1. Based on: run52
;; 2. Description: 8 ; M3 with covariates less IIV Dummy

;; x1. Author: user

$PROBLEM    43 ;1 compartment pharmacokinetic model for artesunate and dihydroartemisinin
$INPUT      ID ;ID
            DATE=DROP ;Date if available
            TIME ;Time if date is present this chould be clock time otherwise time in hours
            ODV ;Dependent variable (nonlogarithmic in this example) 
            DV ;Dependent variable (logarithmic in this example)
            OODV ;Dependent variable (logarithmic in this example,including values below limit of quantification) 
            WT ;Body weight (covariate) 
            EVID ;Event ID record 
            MDV1 ;Missing dependent variable (1=dose and concentrations below limit of quantification)
            MDV ;Missing dependent variable (1=dose)
            AMT ;Dose amount
            FLAG ;Flag for artesunate=1,dihydroartemisinin=2
            CMT ;Compartment (1=dose,2=artesunate,3=dihydroartemisinin)
            CMT2 ;Compartment (1=dose and artesunate,2=dihydroartemisinin)  
            BQL ;Below quantification limit (1=BQL) (optional)
            PREG ;Pregnancy (covariate 1=pregnant) 
            PARA ;Parasite count (covariate) 
            LNPC ;Parasite count (covariate,in logarithmic scale)
            SORT ;Sort row
            HB AST ALT BIL EGA ;months
$DATA      Simulated_dataset.csv IGNORE=# IGNORE(EVID.EQ.2)
 ;DOSE(NMOL)CP(NMOL/L)
$SUBROUTINE ADVAN5 TRANS1
$MODEL      COMP=(1) ;(DEPOT DOSE)
            COMP=(2) ;(CONCENTRATION ARTESUNATE)
            COMP=(3) ;(CONCENTRATION DIHYDROARTEMISININ (METABOLITE))
            COMP=(4) ;(TRANSIT COMPARTMENT 1)
            COMP=(5) ;(TRANSIT COMPARTMENT 2)
            COMP=(6) ;(TRANSIT COMPARTMENT 3)
$PK

;;; F1LNPC-DEFINITION START
F1LNPC = ( 1 + THETA(9)*(LNPC - 5.88))
;;; F1LNPC-DEFINITION END


;;; F1ALT-DEFINITION START
F1ALT = ( 1 + THETA(8)*(ALT - 20.75))
;;; F1ALT-DEFINITION END

;;; F1-RELATION START
F1COV=F1ALT*F1LNPC
;;; F1-RELATION END


;;; CLMPREG-DEFINITION START
IF(PREG.EQ.1) CLMPREG = 1  ; Most common
IF(PREG.EQ.0) CLMPREG = ( 1 + THETA(7))
;;; CLMPREG-DEFINITION END

;;; CLM-RELATION START
CLMCOV=CLMPREG
;;; CLM-RELATION END


TVCLP  = THETA(1)*((WT/52)**0.75)       ; Population artesunate clearance, allometrically scaled by median body weight with a factor of 0.75
CLP    = TVCLP*EXP(ETA(1))              ; Individual estimate of artesunate clearance 

TVV2   = THETA(2)*((WT/52)**1)          ; Population artesunate volume of distribution, allometrically scaled by median body weight with a factor of 1
V2     = TVV2*EXP(ETA(2))               ; Individual estimate of artesunate volume of distribution

TVCLM  = THETA(3)*((WT/52)**0.75)   ; Population dihydroartemisinin clearance, allometrically scaled by median body weight with a factor of 0.75

TVCLM = CLMCOV*TVCLM
CLM    = TVCLM*EXP(ETA(3))              ;

TVV3   = THETA(4)*((WT/52)**1)          ; Population dihydroartemisinin volume of distribution, allometrically scaled by median body weight with a factor of 1
V3     = TVV3*EXP(ETA(4))               ;

TVMT   = THETA(5)                       ; Population artesunate mean transit time
MT     = TVMT*EXP(ETA(5))               ; Individual estimate of artesunate mean transit time

TVF1   = THETA(6)                       ; 

TVF1 = F1COV*TVF1
F1     = TVF1*EXP(ETA(6))               ; Artesunate relative bioavailability

S2=V2                                   ; Scaling factor for artesunate
S3=V3                                   ; Scaling factor for dihydroartemisinin

NN    = 3 ;Number of transit compartments
KTR   = (NN+1)/MT ;
K14   = KTR
K45   = KTR
K56   = KTR
K62   = KTR

K23=CLP/V2
K30=CLM/V3

;IF(CMT.EQ.2) LLOQ=1.1383         ;value of LLOQ for artesunate in nmol/L (log)
;IF(CMT.EQ.3) LLOQ=1.9509         ;value of LLOQ for dihydroartemisinin nmol/L (log)

HTA = LOG(2)/K23 ;Half-life artesunate
HTD = LOG(2)/K30 ;Half-life dihydroartemisinin

$ERROR
                                                                  ;M3-method for evaluating values below LLOQ
  IF(CMT.EQ.2) THEN
    IPRED    = A(2)/V2                                                  ;Predicted plasma concentration of artesunate
    W     = SQRT(SIGMA(1,1))                                            ;Additive residual error on log scale
  ENDIF

  IF(CMT.EQ.3) THEN
    IPRED    = A(3)/V3                                                  ;Predicted plasma concentration of dihydroartemisinin
    W     = SQRT(SIGMA(2,2))                                            ;Additive residual error on log scale
  ENDIF

  IF(IPRED.GT.0)  IPRED = LOG(IPRED)
 IF (CMT.EQ.2) LLOQ = 1.1383                                           ;Value of artesunate LLOQ
  IF (CMT.EQ.3) LLOQ = 1.9509                                           ;Value of dihydroartemisinin LLOQ

  DUM  = (LLOQ-IPRED)/W                                                 ;Positive when IPRED is larger then LLOQ
  CUMD = PHI(DUM)

;-- Prediction DV>=LOQ ------------------------------------------------------

  IRES  = IPRED-DV
  IWRES = IRES/W

  IF(BQL.EQ.0.AND.CMT.EQ.2) THEN
      F_FLAG = 0
      Y      = IPRED+ERR(1)
  ENDIF

  IF(BQL.EQ.0.AND.CMT.EQ.3) THEN
      F_FLAG = 0
      Y      = IPRED+ERR(2)
  ENDIF

;-- Likelihood DV<LOQ -------------------------------------------------------

   IF(BQL.EQ.1) THEN
      F_FLAG = 1
      Y      = CUMD + 0.000001
   ENDIF
  IF(AMT.GT.0)    DTIM  = TIME
                  TAD   = TIME-DTIM

$THETA
(0, 3570) ; 1.CLP
(0, 1700) ; 2.V2
(0, 190) ; 3.CLM
(0, 267) ; 4.V3
(0, 0.832) ; 5.MTT
(1) FIX ; 6.F1
(-1, -0.214) ; CLMPREG1
(-0.024, 0.0215,0.057) ; F1ALT1
(-0.199, 0.138,0.334) ; F1LNPC1

$OMEGA  
 0.0672 ;       1.CL
 0 FIX  ;      2.V2_
 0.00809 ;     3.CLM_
 0 FIX  ;       4.VM
 0.32 ;     5.MTT_
 0.0887 ;       8.F1

$SIGMA 0.892 ;    RUV_ARS
 0.66 ;    RUV_DHA

$EST
MAXEVAL=0 PRINT=5 POSTHOC METHOD=1 LAPLACIAN INTER
;$COV
$TABLE ID TIME TAD CLP V2 CLM V3 MT F1 HTA HTD CMT FLAG IPRED
PRED CWRES MDV NOPRINT ONEHEADER
FILE=mytab1

