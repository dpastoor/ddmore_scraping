;Jönsson S, Davidse A, Wilkins J, Van der Walt JS, Simonsson US, Karlsson MO, Smith P, McIlleron H.
;Population pharmacokinetics of ethambutol in South African tuberculosis patients.
;Antimicrob Agents Chemother. 2011 Sep;55(9):4230-7. doi: 10.1128/AAC.00274-11.
$PROBLEM EMB (log) 
$INPUT ID STID=DROP DAY=DROP RATE=DROP OT1=DROP OT2=DROP TIME ODV=DROP DV MDV AMT
       ADDL II EVID AGE=DROP SEX=DROP WT HT=DROP BMI=DROP RACE=DROP SMOK=DROP ALC=DROP
       DPKG=DROP HIV HB=DROP HCT=DROP RBC=DROP MCV=DROP WBC=DROP AP=DROP LNAP=DROP
       ALT=DROP AST=DROP CRT=DROP TBIL=DROP UREA=DROP FDC=DROP LOC=DROP CLCR=DROP BSA=DROP
       INHP=DROP PZAP=DROP CLC2=DROP OCC DOS=DROP NDOS=DROP DTIM=DROP TAD=DROP
       TAD2=DROP
$DATA Simulated_data.csv IGNORE=@ 
;TIME in hours
;DV is logtransformed concentrations
;MDV Missing dependent variable: 0 for observations, 1 for dose records
;AMT dose amount
;ADDL additional dose data items, defines number of doses given in addition to the current
;II interdose interval, the dosing interval between repeated doses
;EVID event data item: 0 for observations, 1 for dose records
;WT body weight in kg
;HIV: 0 HIV negative, 1 HIV positive

$SUBROUTINE ADVAN5 TRANS1

$MODEL NCOMP=4
  COMP = (TRANSIT)
  COMP = (ABS)
  COMP = (CENTRAL)
  COMP = (PERIPH)

$PK
    OC1   = 0
    OC2   = 0
    OC3   = 0
    OC4   = 0

    IF(OCC.EQ.1) OC1 = 1
    IF(OCC.EQ.2) OC2 = 1
    IF(OCC.EQ.3) OC3 = 1
    IF(OCC.EQ.4) OC4 = 1

    TVCL  = THETA(1)*(WT/50)**0.75     ;allometric scaling with body weight,
    TVV2  = THETA(2)*(WT/50)**1        ;on clearance and volume terms. Theory based exponents
    TVKA  = THETA(3)
    TVV3  = THETA(6)*(WT/50)**1
    TVQ   = THETA(7)*(WT/50)**0.75
    TVMTT = THETA(8)

    FCOV = 1
    IF(HIV.EQ.1) FCOV = 1 + THETA(9)   ;HIV on bioavailability

    CL    = TVCL*EXP(ETA(1)+OC1*ETA(7)+OC2*ETA(8)+OC3*ETA(9)+OC4*ETA(10))
    V2    = TVV2*EXP(ETA(2))
    KA    = TVKA*EXP(ETA(3))
    V3    = TVV3*EXP(ETA(4))
    Q     = TVQ*EXP(ETA(5))
    F1    = 1*FCOV
    MTT = TVMTT * EXP(ETA(6))
    KTR = 1/MTT

  K12 = KTR
  K23 = KA
  K34 = Q/V2
  K43 = Q/V3
  K30 = CL/V2

  S3 = V2

$ERROR (ONLY OBSERVATIONS)     ;log transformation both sides
    IF (F.EQ.0) THEN
      IPRED = 0
    ELSE
      IPRED = LOG(F)
    ENDIF
    IRES = DV - IPRED
    W = SQRT(THETA(4)**2 + (THETA(5)/(F+.0001))**2)
    IF (W.EQ.0) W = 1
    IWRES = IRES / W
    Y = IPRED + W * EPS(1)

$THETA  (0,40.90000) ; 1CL
 (0,139.0000)        ; 2V
 (0,0.505000)        ; 3KA
 (0,0.314000)        ; 4Prop error
 (0,0.120000)        ; 5Add error
 (0,1110.000)        ; 6V3
 (0,28.50000)        ; 7Q
 (0,0.359000)        ; 8MTT
 (-1,-0.15500)       ; 9FCOV
$OMEGA  0.040800     ; 1CL
 0  FIX              ; 2V2
 0.366000            ; 3KA
 0  FIX              ; 4V3
 0  FIX              ; 5Q
 0.702000            ; 6MTT
$OMEGA  BLOCK(1)
 0.203000            ; 7IOV-CL
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME

$SIGMA  1  FIX
$ESTIMATION METHOD=1 INTER NOABORT PRINT=2 MAXEVAL=0
;$COVARIANCE



