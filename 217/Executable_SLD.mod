; Claret-like model for describing CTS. No resistance to drug treatment. Independent additive effect of the two drugs. M3 method for BLQ. Additive error

$SIZES LIM6=2000
$PROBLEM Model for SLD(t)

$INPUT CCOM,ID,
TIME, ;in day
DV, ;SLD in mm
BQL, ;the LLOQ is 5 mm
CB, ;exposure to carboplatin (per-cycle average AUC)
G, ;exposure to gemcitabine (per-cycle average AUC)
EVID,
FLG, ;FLG=2 for SLD data, FLG=1 for exposure-related entries
CMT ;CMT=1 for SLD data

$DATA Simulated_SLD.csv IGNORE=C 

$SUBROUTINE ADVAN6 TOL=3

$MODEL
;Tumour
COMP = (TUMOUR, DEFOBS)

$PK
KG = THETA(1)*EXP(ETA(1)) ;tumour growth rate constant (1/day)
KD0	= THETA(2)*EXP(ETA(2)) ;carboplatin related death rate constant (1/day/AUC0)
KD1 = THETA(3)*EXP(ETA(2)) ;gemcitabine related death rate constant (1/day/AUC1)
IBASE = THETA(4)*EXP(ETA(3)) ;baseline SLD (m)
FADD = THETA(5) ;SD of additive error (mm)

; ==== SLD baseline ====
A_0(1) = IBASE*1000

; Backward interpolation of exposure-related data
IF(NEWIND.NE.2) OCB=CB
IF(NEWIND.NE.2) OG=G
E0 = OCB 
E1 = OG 
OCB=CB
OG=G

$DES
; Model for dSLD(t) 
DADT(1) = KG/1000 * A(1) - (KD0/1000 * E0 + KD1/100 * E1) * A(1) 

$ERROR
LLOQ = 5 ;5 mm is LLOQ
IPRED = A(1)

W = FADD ;SD of additive unexpained variability

; Probability of SLD<LLOQ
DUM = (LLOQ-IPRED)/W
DUM2 = PHI(DUM) 

IF(BQL.EQ.1) THEN
	F_FLAG = 1
	Y = DUM2
ENDIF
IF(BQL.EQ.0) THEN
	F_FLAG = 0
	Y = IPRED+ERR(1)*W
ENDIF

IRES = IPRED-DV
IWRES = IRES/W

;Parameters
$THETA  
(0, 0.3) ; KG [1/day]
(0, 0.03) ; KD0 [1/day/AUC0]
(0,0.01) ; KD1 [1/day/AUC1]
(0, 0.065) ;IBASE [m]
(0, 20) ; FADD 

$OMEGA 
0.08 ; KG
0.1 ; KD
0.1 ;IBASE

;Error
$SIGMA  1 FIX  ; placeholder

$EST MAXEVAL=9000 PRINT=10 METH=1 LAPLACIAN INTER NUMERICAL SLOW NOABORT
$COV SLOW

$TABLE ID TIME DV IPRED PRED IRES IWRES W CWRES EVID ETA1 ETA2 ETA3 KG KD0 KD1 IBASE 
FILE=sdtab_SLD NOPRINT ONEHEADER
