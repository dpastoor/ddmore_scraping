$SIZES      NO=1000 LIM6=1000 LNP4=16000 DIMTMP=6000 LVR=50
$PROBLEM    OS - IL_SUV - SLD
$INPUT      
 ID    ; Individual patient identifier
 TIME  ; Time in hours
 FLAG  ; Record identifier: 1: dosing record, 3=overall survival; 4=SLD; 13=censoring from OS; 15=individual lesion SUVmax
 DV    ; dependent variable (log-transformed for individual lesion SUVmax)
 CL    ; post-hoc clearance from previously-developed PK model
 EVID  ; event identifier
 CENS  ; overall survival/censoring information: 0=censored; 1=event
 DOS   ; daily dose in mg
 LESI  ; lesion number for SUV
 L2    ; L2 data item

$DATA       Simulated_SLD_SUV_OS_GIST.csv IGNORE=#

$SUBROUTINE ADVAN6 TOL=6
$MODEL      NCOMPS=9

$PK
; --- Sunitinib daily AUC
AUC     =    DOS/CL

; ---- Initialize flags and variables used in $ERROR
IF(NEWIND.NE.2) THEN
; -- Flag for presence or absence of lesion 1-5 for SUVmax
ISLES1 = 0
ISLES2 = 0
ISLES3 = 0
ISLES4 = 0
ISLES5 = 0
; -- Amount in compartment 1-5 at week 1 (168h) (SUVmax lesion 1-5)
WK1A1 = 0
WK1A2 = 0
WK1A3 = 0
WK1A4 = 0
WK1A5 = 0
; -- Relative change from baseline at week 1 in compartments 1-5 (SUVmax lesion 1-5)
RCFB1A1 = 0
RCFB1A2 = 0
RCFB1A3 = 0
RCFB1A4 = 0
RCFB1A5 = 0
; -- Maximum relative change in SUV from baseline at week 1 across lesions
RCFB1MAX= 0
ENDIF

;=======================================
; Parameters and variables for SUVmax model
;=======================================
;----- Parameters for IDR model --------------
; ~~~~~~~~ Kout: rate constant for the loss of response
; (here all OMEGAs fixed to 0 as no IIV or ILV was identified)
TVKOUT  =   THETA(1)/1000/24/7
KOUT1   =   TVKOUT*EXP(ETA(1)+ETA(2))      ; lesion 1
KOUT2   =   TVKOUT*EXP(ETA(1)+ETA(3))      ; lesion 2
KOUT3   =   TVKOUT*EXP(ETA(1)+ETA(4))      ; lesion 3
KOUT4   =   TVKOUT*EXP(ETA(1)+ETA(5))      ; lesion 4
KOUT5   =   TVKOUT*EXP(ETA(1)+ETA(6))      ; lesion 5

; ~~~~~~~~ Slope of disease progression 
; (here THETA and all OMEGAs fixed to 0 as no disease progression was identified)
TVALPHA =   THETA(2)/1000/24/7
ALPHA1  =   TVALPHA*EXP(ETA(7)+ETA(8))     ; lesion 1
ALPHA2  =   TVALPHA*EXP(ETA(7)+ETA(9))     ; lesion 2
ALPHA3  =   TVALPHA*EXP(ETA(7)+ETA(10))    ; lesion 3
ALPHA4  =   TVALPHA*EXP(ETA(7)+ETA(11))    ; lesion 4
ALPHA5  =   TVALPHA*EXP(ETA(7)+ETA(12))    ; lesion 5

; ~~~~~~~~ Sunitinib drug effect parameter on SUVmax 
;(DRUG_SUV in original publication)
TVDRUG  =   THETA(3)
DRUG1   =   TVDRUG*EXP(ETA(14)+ETA(15))    ; lesion 1
DRUG2   =   TVDRUG*EXP(ETA(14)+ETA(16))    ; lesion 2
DRUG3   =   TVDRUG*EXP(ETA(14)+ETA(17))    ; lesion 3
DRUG4   =   TVDRUG*EXP(ETA(14)+ETA(18))    ; lesion 4
DRUG5   =   TVDRUG*EXP(ETA(14)+ETA(19))    ; lesion 5

; ~~~~~~~~ Baseline SUVmax 
; (SUV_0 in original publication)
TVBASE  =   THETA(4)
IBASE1  =   TVBASE*EXP(ETA(21)+ETA(22))   ; lesion 1
IBASE2  =   TVBASE*EXP(ETA(21)+ETA(23))   ; lesion 2
IBASE3  =   TVBASE*EXP(ETA(21)+ETA(24))   ; lesion 3
IBASE4  =   TVBASE*EXP(ETA(21)+ETA(25))   ; lesion 4
IBASE5  =   TVBASE*EXP(ETA(21)+ETA(26))   ; lesion 5

;----- Parameters for Cholesky decomposition --
; ~~~~~~~~ Diagonal (Variance)
VD11    =   THETA(5)
VD22    =   THETA(5)
VD33    =   THETA(5)
VD44    =   THETA(5)
VD55    =   THETA(5)

; ~~~~~~~~ Off-diagonal (Correlation)
CR12    =   THETA(6)
CR13    =   THETA(6)
CR14    =   THETA(6)
CR15    =   THETA(6)
CR23    =   THETA(6)
CR24    =   THETA(6)
CR25    =   THETA(6)
CR34    =   THETA(6)
CR35    =   THETA(6)
CR45    =   THETA(6)

;=======================================
; Parameters and variables for SLD model
;=======================================
TVKG    =    THETA(7)/24/7/1000
KG      =    TVKG  *EXP(ETA(27))     ; Tumor growth rate constant (K_GROW)
TVLAM   =    THETA(8)/24/7/1000
LAMBDA  =    TVLAM *EXP(ETA(28))     ; Tumor size reduction rate constant
TVDRU   =    THETA(9)/24/7/1000
KDRUG   =    TVDRU *EXP(ETA(13))     ; Tumor size reduction rate constant (DRUG_SLD)
IBASE   =    THETA(10)*EXP(ETA(20))  ; Baseline tumor size (SLD_0)
KEO     =    LOG(2)/50               ; Equilibration rate constant for the effect compartment

;=======================================
; Parameters for time-to-event models
;=======================================
LAMBH     =   THETA(12)/24/7         ; Scale parameter in the Weibull pdf for the survival model
ALPHH     =   THETA(13)              ; Shape parameter in the Weibull pdf for the survival model
LAMBD     =   THETA(14)/24/7         ; Scale parameter in the Weibull pdf for the drop out model
ALPHD     =   THETA(15)              ; Shape parameter in the Weibull pdf for the drop out model

;-----Compartment initialization-----
A_0(1)  =   IBASE1                   ; SUVmax lesion 1
A_0(2)  =   IBASE2                   ; SUVmax lesion 2
A_0(3)  =   IBASE3                   ; SUVmax lesion 3
A_0(4)  =   IBASE4                   ; SUVmax lesion 4
A_0(5)  =   IBASE5                   ; SUVmax lesion 5

A_0(6)  =   IBASE                    ; SLD

$DES

;=======================================
;             SUV model
;=======================================
; ~~~~~~~~ Kin: rate constant for the production of SUVmax response
KIN1    =   KOUT1*IBASE1                   ; lesion 1
KIN2    =   KOUT2*IBASE2                   ; lesion 2
KIN3    =   KOUT3*IBASE3                   ; lesion 3
KIN4    =   KOUT4*IBASE4                   ; lesion 4
KIN5    =   KOUT5*IBASE5                   ; lesion 5

; ~~~~~~~~ Drug effect
EFF1    =   DRUG1*A(9)                     ; lesion 1
EFF2    =   DRUG2*A(9)                     ; lesion 2
EFF3    =   DRUG3*A(9)                     ; lesion 3
EFF4    =   DRUG4*A(9)                     ; lesion 4
EFF5    =   DRUG5*A(9)                     ; lesion 5

; ~~~~~~~~ Differential equations
DADT(1) =   KIN1*(1+ALPHA1*T)-KOUT1*(1+EFF1)*A(1)      ; lesion 1
DADT(2) =   KIN2*(1+ALPHA2*T)-KOUT2*(1+EFF2)*A(2)      ; lesion 2
DADT(3) =   KIN3*(1+ALPHA3*T)-KOUT3*(1+EFF3)*A(3)      ; lesion 3
DADT(4) =   KIN4*(1+ALPHA4*T)-KOUT4*(1+EFF4)*A(4)      ; lesion 4
DADT(5) =   KIN5*(1+ALPHA5*T)-KOUT5*(1+EFF5)*A(5)      ; lesion 5
DADT(9) =   KEO*(AUC-A(9))                             ; Effect compartment

;=======================================
;             SLD model
;=======================================
AUC1    =   A(9)*KDRUG                                 ; exposure driven drug effect
DADT(6) =   KG*A(6)- AUC1*EXP(-(LAMBDA*T))*A(6)        ; longitudinal model describing tumor growth

;=======================================
;             Time-to-event model
;=======================================
; ~~~~~~~~ Weibull model for overall survival
DELX=1E-16
DADT(7) = LAMBH*ALPHH*(T+DELX)**(ALPHH-1)* EXP(THETA(16)*RCFB1MAX)
; ALPHH fixed to 1 to get constant hazard instead of Weibull

; ~~~~~~~~ Weibull model for censoring
DADT(8) = LAMBD*ALPHD*(T+DELX)**(ALPHD-1)
; ALPHH fixed to 1 to get constant hazard instead of Weibull


$ERROR  
;=======================================
;              SUVmax model
;=======================================
;---- Cholesky decomposition ----------------
; ~~~~~~~~ Off-diagonal (Covariance)
CV12    =   SQRT(VD11*VD22)*CR12
CV13    =   SQRT(VD11*VD33)*CR13
CV14    =   SQRT(VD11*VD44)*CR14
CV15    =   SQRT(VD11*VD55)*CR15
CV23    =   SQRT(VD22*VD33)*CR23
CV24    =   SQRT(VD22*VD44)*CR24
CV25    =   SQRT(VD22*VD55)*CR25
CV34    =   SQRT(VD33*VD44)*CR34
CV35    =   SQRT(VD33*VD55)*CR35
CV45    =   SQRT(VD44*VD55)*CR45

; ~~~~~~~~ Cholesky decomposition matrix
CD11    =   SQRT(VD11)
CD12    =   CV12 / CD11
CD13    =   CV13 / CD11
CD14    =   CV14 / CD11
CD15    =   CV15 / CD11

CD22    =   SQRT(VD22 - (CD12*CD12))
CD23    =   (CV23 - CD13*CD12) / CD22
CD24    =   (CV24 - CD14*CD12) / CD22
CD25    =   (CV25 - CD15*CD12) / CD22

CD33    =   SQRT(VD33 - (CD13*CD13 + CD23*CD23))
CD34    =   (CV34 - (CD14*CD13 + CD23*CD24)) / CD33
CD35    =   (CV35 - (CD15*CD13 + CD23*CD25)) / CD33

CD44    =   SQRT(VD44 - (CD14*CD14 + CD24*CD24 + CD34*CD34))
CD45    =   (CV45 - (CD14*CD15 + CD24*CD25 + CD34*CD35)) / CD44

CD55    =   SQRT(VD55 - (CD15*CD15 + CD25*CD25 + CD35*CD35 + CD45*CD45))

EPS1    =   CD11*EPS(1)
EPS2    =   CD12*EPS(1)+CD22*EPS(2)
EPS3    =   CD13*EPS(1)+CD23*EPS(2)+CD33*EPS(3)
EPS4    =   CD14*EPS(1)+CD24*EPS(2)+CD34*EPS(3)+CD44*EPS(4)
EPS5    =   CD15*EPS(1)+CD25*EPS(2)+CD35*EPS(3)+CD45*EPS(4)+CD55*EPS(5)

; --- Individual lesion SUVmax
IF(FLAG.EQ.15) IPRED=0

IF(FLAG.EQ.15.AND.LESI.EQ.1.AND.A(1).GT.0) IPRED = LOG(A(1))                 ; lesion 1
IF(FLAG.EQ.15.AND.LESI.EQ.2.AND.A(2).GT.0) IPRED = LOG(A(2))                 ; lesion 2
IF(FLAG.EQ.15.AND.LESI.EQ.3.AND.A(3).GT.0) IPRED = LOG(A(3))                 ; lesion 3
IF(FLAG.EQ.15.AND.LESI.EQ.4.AND.A(4).GT.0) IPRED = LOG(A(4))                 ; lesion 4
IF(FLAG.EQ.15.AND.LESI.EQ.5.AND.A(5).GT.0) IPRED = LOG(A(5))                 ; lesion 5

IF(FLAG.EQ.15.AND.LESI.EQ.1) EPSN = EPS1                                     ; lesion 1
IF(FLAG.EQ.15.AND.LESI.EQ.2) EPSN = EPS2                                     ; lesion 2
IF(FLAG.EQ.15.AND.LESI.EQ.3) EPSN = EPS3                                     ; lesion 3
IF(FLAG.EQ.15.AND.LESI.EQ.4) EPSN = EPS4                                     ; lesion 4
IF(FLAG.EQ.15.AND.LESI.EQ.5) EPSN = EPS5                                     ; lesion 5

IF(FLAG.EQ.15) THEN
  IRES  = DV-IPRED
  W     = 1
  IWRES = IRES/W
  F_FLAG = 0                       ; prediction
  Y = IPRED+EPSN
ENDIF

;=======================================
;              SLD model
;=======================================
IF(FLAG.EQ.4) THEN                 ; SLD
  IPRED   =   A(6)
  IRES    =   DV-IPRED
  W1      =   IPRED*THETA(11)
ENDIF
IF (W1.EQ.0) W1 = 1
IF(FLAG.EQ.4) THEN
  IWRES   =   IRES/W1
  F_FLAG  =   0                    ; prediction
  Y       =   IPRED+W1*EPS(6)
ENDIF

;=======================================
;             Time-to-event model
;=======================================
; --- Flag for presence of lesion in dataset; switches from 0 to 1 if lesion present
IF(LESI.EQ.1) ISLES1 = 1
IF(LESI.EQ.2) ISLES2 = 1
IF(LESI.EQ.3) ISLES3 = 1
IF(LESI.EQ.4) ISLES4 = 1
IF(LESI.EQ.5) ISLES5 = 1

; --- Individual SUVmax prediction conditional on presence of lesion
AA1 = A(1)*ISLES1
AA2 = A(2)*ISLES2
AA3 = A(3)*ISLES3
AA4 = A(4)*ISLES4
AA5 = A(5)*ISLES5

; --- Lesion with largest relative change in SUV from baseline at week 1
IF(FLAG.EQ.1.AND.TIME.EQ.168) THEN
 WK1A1 = AA1                                ; SUV at week 1 for lesion 1
 WK1A2 = AA2                                ; SUV at week 1 for lesion 2
 WK1A3 = AA3                                ; SUV at week 1 for lesion 3
 WK1A4 = AA4                                ; SUV at week 1 for lesion 4
 WK1A5 = AA5                                ; SUV at week 1 for lesion 5
 
 RCFB1A1 = (WK1A1-IBASE1)/IBASE1*ISLES1     ; relative change in SUV from baseline at week 1 - lesion 1
 RCFB1A2 = (WK1A2-IBASE2)/IBASE2*ISLES2     ; relative change in SUV from baseline at week 1 - lesion 2
 RCFB1A3 = (WK1A3-IBASE3)/IBASE3*ISLES3     ; relative change in SUV from baseline at week 1 - lesion 3
 RCFB1A4 = (WK1A4-IBASE4)/IBASE4*ISLES4     ; relative change in SUV from baseline at week 1 - lesion 4
 RCFB1A5 = (WK1A5-IBASE5)/IBASE5*ISLES5     ; relative change in SUV from baseline at week 1 - lesion 5
ENDIF

; maximum relative change in SUVmax from baseline across lesions, at week 1
IF(FLAG.EQ.1.AND.TIME.EQ.168.AND.RCFB1A1.LT.RCFB1A2.AND.RCFB1A1.LT.RCFB1A3.AND.RCFB1A1.LT.RCFB1A4.AND.RCFB1A1.LT.RCFB1A5) RCFB1MAX = RCFB1A1
IF(FLAG.EQ.1.AND.TIME.EQ.168.AND.RCFB1A2.LT.RCFB1A1.AND.RCFB1A2.LT.RCFB1A3.AND.RCFB1A2.LT.RCFB1A4.AND.RCFB1A2.LT.RCFB1A5) RCFB1MAX = RCFB1A2
IF(FLAG.EQ.1.AND.TIME.EQ.168.AND.RCFB1A3.LT.RCFB1A1.AND.RCFB1A3.LT.RCFB1A2.AND.RCFB1A3.LT.RCFB1A4.AND.RCFB1A3.LT.RCFB1A5) RCFB1MAX = RCFB1A3
IF(FLAG.EQ.1.AND.TIME.EQ.168.AND.RCFB1A4.LT.RCFB1A1.AND.RCFB1A4.LT.RCFB1A2.AND.RCFB1A4.LT.RCFB1A3.AND.RCFB1A4.LT.RCFB1A5) RCFB1MAX = RCFB1A4
IF(FLAG.EQ.1.AND.TIME.EQ.168.AND.RCFB1A5.LT.RCFB1A1.AND.RCFB1A5.LT.RCFB1A2.AND.RCFB1A5.LT.RCFB1A3.AND.RCFB1A5.LT.RCFB1A4) RCFB1MAX = RCFB1A5


; --- TTD (Time to death)--------------
CHZ     =   A(7)           ; cumulative hazard for survival
DCHZ    =   A(8)           ; cumulative hazard for drop out

SUR     =   EXP(-CHZ)      ; survival probability
SURD    =   EXP(-DCHZ)     ; drop out probability

DEL     =   1E-16
HAZNOW  =   0
HAZDNOW =   0
HAZNOW  = LAMBH*ALPHH*(TIME+DEL)**(ALPHH-1)* EXP(THETA(16)*RCFB1MAX)
HAZDNOW = LAMBD*ALPHD*(TIME+DEL)**(ALPHD-1)

IF(FLAG.EQ.3.AND.EVID.EQ.0.AND.CENS.EQ.0) THEN    ; survival
F_FLAG = 1                                      ; prediction
Y = SUR                  ; probability of survival (censored event)
ENDIF

IF(FLAG.EQ.3.AND.EVID.EQ.0.AND.CENS.EQ.1) THEN    ; survival
F_FLAG = 1                                      ; prediction
Y = SUR*HAZNOW           ; probability of event (death) at time=TIME
ENDIF

IF(FLAG.EQ.13.AND.CENS.EQ.1) THEN                 ; censored
F_FLAG = 1                                      ; likelihood
Y = SURD                 ; probability of not dropping out
ENDIF

IF(FLAG.EQ.13.AND.CENS.EQ.0) THEN                 ; censored
F_FLAG = 1                                      ; likelihood
Y = SURD*HAZDNOW         ; probability of dropping out at time=TIME
ENDIF

;-----INITIAL ESTIMATES THETAs -----------------
; ----------- SUV model  --------------
$THETA  555.634 FIX   ; 1  KOUT
$THETA  0 FIX         ; 2  ALPHA
$THETA  0.94569 FIX   ; 3  DRUG
$THETA  7.58866 FIX   ; 4  IBASE_SUV
$THETA  0.173794 FIX  ; 5  VAR RUV
$THETA  0.466827 FIX  ; 6  CORR RUV
; ----------- SLD model  --------------
$THETA  10.4796 FIX   ; 7  KG*1000
$THETA  20.0529 FIX   ; 8  LAMBDA*1000
$THETA  16.5971 FIX   ; 9  KDRUG*1000
$THETA  263.064 FIX   ; 10 IBASE_SLD
$THETA  0.0667619 FIX ; 11 RUV

; ----------- TTE model  --------------
$THETA  (0,0.018792)  ; 12  LAMBH
$THETA  1 FIX         ; 13  ALPHH
$THETA  (0,0.02)      ; 14  LAMBD
$THETA  1 FIX         ; 15  ALPHD
$THETA  5.17          ; 16  Predictor

;-----INITIAL ESTIMATES OMEGAs -----------------
$OMEGA  0  FIX  ; 1  KOUT IIV
$OMEGA  0  FIX  ; 2  KOUT LES1
$OMEGA  0  FIX  ; 3  KOUT LES2
$OMEGA  0  FIX  ; 4  KOUT LES3
$OMEGA  0  FIX  ; 5  KOUT LES4
$OMEGA  0  FIX  ; 6  KOUT LES5

$OMEGA  0  FIX  ; 7  ALPHA IIV
$OMEGA  0  FIX  ; 8  ALPHA LES1
$OMEGA  0  FIX  ; 9  ALPHA LES2
$OMEGA  0  FIX  ; 10 ALPHA LES3
$OMEGA  0  FIX  ; 11 ALPHA LES4
$OMEGA  0  FIX  ; 12 ALPHA LES5

$OMEGA  BLOCK(2)  FIX
 0.398557             ; 13 KDRUG SLD
 0.397941 0.548112    ; 14 DRUG IIV

$OMEGA  BLOCK(1)  FIX
 0.324481             ; 15 DRUG LES1
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME

$OMEGA  0.288669  FIX ; 20 IBASE SLD
$OMEGA  0.105177  FIX ; 21 IBASE IIV

$OMEGA  BLOCK(1)  FIX
 0.0549036            ; 22 IBASE LES1
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME

$OMEGA  0.361258  FIX ; 27 KG
$OMEGA  0  FIX        ;  28 LAMBDA

$SIGMA
 1  FIX  ;   lesion 1
 1  FIX  ;   lesion 2
 1  FIX  ;   lesion 3
 1  FIX  ;   lesion 4
 1  FIX  ;   lesion 5
 1  FIX  ;   SLD

$ESTIMATION MAXEVAL=0 METHOD=1 LAPLACIAN NOABORT PRINT=1 NSIG=2 SIGL=6
;$COV UNCONDITIONAL MATRIX=R PRINT=E
