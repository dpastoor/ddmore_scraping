PROB NIVOLUMAB Population PK Analysis
; DESCRIPTION: Linear 2-cmpt PK model, with time-varying CL (sigmoid-Emax) - Final PPK model
; NOTES: Method: FOCEI; proportional residual error model
; Definitions of Variables in PPK dataset
; TIME: Actual time after first dose
; BBWT: Baseline body weight
; BALB: Baseline Albumin 
; BLDH: Baseline lactate dehydrogenase 
; TTYPEN: Tumor type (1=MEL; 2=NSCLC; 3=OTHERS; 4=RCC)
; HEPAN: Hepatic function (1=normal)
; RACEN: Race (1=White; 2=African American; 3=Asian, 4=Others)
; PS: Performance Status (0, 1)
; SEXN: 1=Male; 2=Female
; FLAGN: Exclusion variable (1=excluded records)

$INPUT  ID TIME EVID AMT RATE DV AGE SEXN BBWT BGFR PS RACEN 
       FLAGN MDV 


$DATA ../pkdata1.csv IGNORE=@ 
 IGNORE=(FLAGN=1)

$SUBROUTINE ADVAN3 TRANS4 

$SIM (NEW 123456) ONLY SIM SUBPROBLEMS=1

$COVARIANCE PRINT=ER

$THETA 
 (0, 0.010) ; CL [L/h] ; theta1
 (0, 4) ; VC [L] ; theta2
 (0, 0.04) ; Q [L/h] ; theta3
 (0, 4) ; VP [L] ; theta4
 ;
 0 FIX ; AERR [ug/mL] ; theta5
 (0, 0.5) ; PERR [-] ; theta6
 ;
(0.05)  ; CL_BBWT ; theta7
(0.1) ; CL_GFR ; theta8
(0.1) ; CL_SEX ; theta9
(0.05) ; CL_PS ; theta10
(-0.05) ; CL_RAAS; theta11
(0.05)  ; VC_BBWT ; theta12
(0.1) ; VC_SEX ; theta13
;
0.01 ; CL_EMAX ; theta14
1800 ; CL_T50 ; theta15
1 ; CL_HILL; theta16

$OMEGA BLOCK(2)
0.5 ; ZCL [-]
0.1 0.3 ; ZVC [-]
;
$OMEGA
0.3 ; ZVP [-]
0.1 ; ZEMAX [h]

$SIGMA  
 1 FIX ; RESERR

$PK 
; Define fixed effects parameters
 ACL = THETA(1)
 AV1 = THETA(2)
 AQ = THETA(3)
 AV2 = THETA(4)

; Define parameters of time-varying CL (sigmoid Emax)
 AEMAX = THETA(14)
 T50 = THETA(15)
 HILL = THETA(16)

; Define parameters of residual error model 
 AERR = THETA(5) ; Additive
 PERR = THETA(6) ; proportional
 
; Define parameters of between-subject variability 
 ZCL = ETA(1)
 ZV1 = ETA(2)
 ZV2 = ETA(3)
 ZEMAX = ETA(4)
 
; Impute covariate values and define effects
; Missing continuous covariates are imputed as median
; Missing categorical covariates are imputed as mode
 
 MVAL = -99 ; missing value in dataset

; Baseline Body Weight
 BBWT_R = 80 ; refernce value, kg
 BBWT_I = BBWT
 IF (BBWT_I .EQ. MVAL) BBWT_I = BBWT_R ; imputation to reference (~median)
;
 CL_BBWT = (BBWT_I/BBWT_R)**THETA(7) 

;Baseline GFR
 BGFR_R = 90 ; reference value, mL/min
 BGFR_I = BGFR
 IF (BGFR_I .EQ. MVAL) BGFR_I = BGFR_R; GFR imputation to median
 CL_BGFR = (BGFR_I/BGFR_R)**THETA(9)

; SEX
 SEX_I = SEXN
 IF (SEX_I .EQ. MVAL) SEX_I = 1; male
 CL_MALE = THETA(9) ; effect of female:male

; Performance Status (PS)
 PS_I = PS
 IF (PS .EQ. MVAL) PS_I = 1 ; impute missing PS to 1 
 CL_PS_1 = THETA(10); effect of PS 1+:0

; Race
 RACE_I = RACEN
 IF (RACE_I .EQ. MVAL) RACE_I = 1; White
 CL_RAAS = THETA(11) ; effect of Asian:White/Other

; VC
 V1_BBWT = (BBWT_I/BBWT_R)**THETA(12)
 V1_SEX = THETA(13)

; Time-varying CL (Sigmoid Emax Function)
 EMAX = AEMAX + ZEMAX
 CL_TIME = EXP(EMAX*TIME**HILL/(T50**HILL+TIME**HILL))

; Define model parameters: add covariate effect
 TVCL = ACL * CL_BBWT  * CL_BGFR  * CL_TIME
 IF (SEX_I .EQ. 1) TVCL = TVCL * EXP(CL_MALE) ; reference is SEX=2 (Female)
 IF (PS_I .GE. 1) TVCL = TVCL * EXP(CL_PS_1) ; reference is PS=0

 IF (RACE_I .EQ. 3) TVCL = TVCL * EXP(CL_RAAS) ; reference is White/Others
 
;  
 TVV1 = AV1 * V1_BBWT
 IF (SEX_I .EQ. 1) TVV1 = TVV1 * EXP(V1_SEX) ; reference is SEX=2 (Female)


; Calculate PK parameters
  CL = TVCL*EXP(ZCL)
  V1 = TVV1*EXP(ZV1)
  Q  = AQ   
  V2 = AV2*EXP(ZV2)
  S1 = V1

; Calculate alpha and beta half-life
 KE = CL/V1
 K23 = Q/V1
 K32 = Q/V2
 AA = KE+K23+K32
 ALPH = (AA+SQRT(AA*AA-4*KE*K32))/2
 BETA = (AA-SQRT(AA*AA-4*KE*K32))/2


 ALPHHL = LOG(2)/ALPH
 BETAHL = LOG(2)/BETA


$ERROR 
; Define Residual Error Weight
 REWT = F*PERR + AERR ; proportional error model (additive error fixed to 0)

 IPRED = F
 IRES = DV - IPRED
 IWRES = IRES/REWT
 Y = IPRED + REWT*EPS(1)


$TABLE ID TIME EVID AMT RATE DV 
       PRED CIPRED RES CIRES CWRES CIWRES  
      CL V1 Q V2 ALPHHL BETAHL EMAX T50 HILL FLAGN 
       NOPRINT ONEHEADER FILE=nm.tab








