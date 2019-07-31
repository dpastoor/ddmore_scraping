;; Pop-PK model for Colistimethate (CMS) and Colistine (Col) during continuous renal replacement therapy
;; Related article : Multicenter population pharmacokinetic study of colistimethate sodium and colistin 
;; dosed as in normal renal function in patients on continuous renal replacement therapy
;; Anne Leuppi-Taegtmeyer, Laurent A Decosterd, Michael Osthoff, Nicolas J Mueller, Thierry Buclin, Natasci Corti
;; Antimicrobial Agents and Chemotherapy 2018
;; Model author: Thierry Buclin

$PROBLEM PK CMS Col in CRRT

$INPUT ID AGE BW SEX EVID QBL QEFF HT ALB DOSE AMT AMTC TINF RATE RATC TIME CMT DV

$DATA Data0.csv IGNORE=#

$SUBROUTINES ADVAN6 TOL=3

$MODEL NCOMP=10 
COMP=(CMSCENTR DEFDOSE)   ; central CMS 
COMP=(COLCENTR)   ; central Col
COMP=(CMSFILTR)   ; filter CMS   
COMP=(COLFILTR)   ; filter Col 
COMP=(CMSCART)    ; cartridge CMS
COMP=(COLCART)    ; cartridge Col
COMP=(CMSEFFL)    ; cumulative effluent CMS
COMP=(COLEFFL)    ; cumulative effluent Col
COMP=(CMSAUC)     ; cumulative AUC of CMS
COMP=(COLAUC)     ; cumulative AUC of Col 

$PK	
TVCL1M = THETA(1) 
CL1M = TVCL1M * EXP(ETA(1))  ; clearance CMS -> Colistin (L/h)
TVV1 = THETA(2) 
V1 = TVV1 * EXP(ETA(2)) ; distribution volume of CMS (L)
K12 = CL1M/V1 
V3 = 0.2   ;  priming volume of filter and lines (L)
IF (HT.EQ.0) THEN
  K13 = (QBL*60/1000)*(1-0.25)/V1
  K31 = (QBL*60/1000)*(1-0.25)/V3
ELSE
  K13 = (QBL*60/1000)*(1-HT)/V1
  K31 = (QBL*60/1000)*(1-HT)/V3
ENDIF
SIV3 = THETA(3) ; sieving of CMS in the absence of albumin value
V5 = 0.3  ; priming volume of cartridge and lines (L)
K35 = SIV3*(QEFF/1000)/V3
K50 = (QEFF/1000)/V5

TVCL2M = THETA(4)
CL2M = TVCL2M * EXP(ETA(3)) ; clearance Colistin -> elimination (L/h)
TVV2 = THETA(5)
V2  = TVV2 * EXP(ETA(4)) ; distribution volume of colistin (L)
K20 = CL2M/V2
V4 = 0.2   ;  priming volume of filter and lines (L)
IF (HT.EQ.0) THEN
  K24 = (QBL*60/1000)*(1-0.25)/V2
  K42 = (QBL*60/1000)*(1-0.25)/V4
ELSE
  K24 = (QBL*60/1000)*(1-HT)/V2
  K42 = (QBL*60/1000)*(1-HT)/V4
ENDIF
SIV4 = THETA(6) ; sieving of CMS
V6 = 0.3  ; priming volume of cartridge and lines (L)
K46 = SIV4*(QEFF/1000)/V4
K60 = (QEFF/1000)/V6

S1 = V1
S2 = V2
S3 = V3
S4 = V4
S5 = V5
S6 = V6
S7 = 1
S8 = 1
S9 = 1
S10 = 1

F1 = 1155.5/1749.8  ; this is to convert AMT (mg of colistimethate Na) into Colistin dose (mg)

; this is to calculate time after LAST dose
IF (AMT.GT.0) THEN
TDOS=TIME
TAD=0.0
ENDIF
IF (AMT.EQ.0) TAD=TIME-TDOS

; differential equations according to model
$DES
DADT(1) = -K12*A(1) -K13*A(1) + K31*A(3)   ; central CMS
DADT(2) = -K24*A(2) -K20*A(2) + K12*A(1) + K42*A(4)   ;  central Col
DADT(3) = K13*A(1) - K31*A(3) - K35*A(3)   ; filter CMS
DADT(4) = K24*A(2) - K42*A(4) - K46*A(4)   ; filter Col
DADT(5) = K35*A(3) - K50*A(5)   ; cartridge CMS
DADT(6) = K46*A(4) - K60*A(6)   ; cartridge Col
DADT(7) = K50*A(5)   ; cumulative effluent CMS
DADT(8) = K60*A(6)   ; cumulative effluent Col
DADT(9) = A(1)/V1    ; cumulative AUC CMS
DADT(10) = A(2)/V2   ; cumulative AUC Col

$ERROR
IPRED = F
    W = SQRT(THETA(7)**2*IPRED**2 + THETA(8)**2)
    Y = IPRED + W*EPS(1)
 IRES = DV-IPRED
IWRES = IRES/W
C1    = A(1)/V1
C2    = A(2)/V2
CFI3    = A(3)/V3
CFI4    = A(4)/V4
CCA5    = A(5)/V5
CCA6    = A(6)/V6
EFF1  = A(7)
EFF2  = A(8)
AUC1  = A(9)
AUC2  = A(10)

$THETA
(0, 2.31)      ; CL1M
(0, 12.1)      ; V1
(0, 1.05,1.2)  ; SIV3
(0, 1.93)      ; CL2M
(0, 70.1)      ; V2
(0, 0.454,1.2) ; SIV4
(0, 0.222)     ; Prop.RE (sd)
(0, 0.459)     ; Add.RE (sd)

$OMEGA
 0.27     ; IIV CL1M
 0.133    ; IIV V1 
 0.304    ; IIV CL2M
 0.247    ; IIV v2

$SIGMA
 1 FIX    ; Residual error

$EST METHOD=1 INTER MAXEVAL=1 NOABORT SIG=3 PRINT=1 POSTHOC  ; just fit the same model
$COV

; Xpose
$TABLE ID TIME TAD EVID DOSE AMT CMT DV MDV PRED RES WRES IPRED IWRES C1 C2 CFI3 CFI4 CCA5 CCA6 EFF1 EFF2 AUC1 AUC2 ONEHEADER NOPRINT NOAPPEND FILE=sdtab00016
$TABLE CL1M V1 CL2M V2 SIV3 SIV4 K12 K13 K31 K20 K24 K42 K35 K46 K50 K60 AGE BW SEX ONEHEADER NOPRINT NOAPPEND FIRSTONLY FILE=patab00016