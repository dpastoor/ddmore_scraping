;; 1. Based on: run1002g
;; 2. Description: Longitudinal Drug Effect P1 - P8; final model;  power 
;; x1. Author: ana


$SIZES      LTH=80 DIMNEW=2000 ;,MAXIDS=1000000
$PROBLEM    Baseline functional system scores
$INPUT STUDY ID DSS=TIME TRT2CD	TRT VISIT EDSS FLAGFS DV RTYPE EVID MDV CREA CRCL CD MCD V1000 DISN NVISIT

$DATA   Simulated_Novakovic_2016_multiplesclerosis_cladribine_irt.csv IGNORE=@

$PRED
; Cholesky decomposition allowing to fix one of the parameters of the OMEGA matrix
; Diagonal (Variance)
OM11 = THETA(3)    ;Random effect 1: Disability Baseline
OM22 = THETA(4)    ;Random efefct 2: slope of progression
OM33 = THETA(5)    ;Random effect 3: Age
OM44 = THETA(6)    ;Random effect 4: MSD
OM55 = THETA(7)    ;Random effect 5: EXNB

; Off-diagonal (Correlation)
COR12 = THETA (8)    ;Dis-slope
COR13 = THETA (9)    ;Dis-Age
COR23 = THETA (10)   ;slope-age
COR14 = THETA (11)   ;Dis-MSD
COR24 = THETA (12)   ;slope-MSD
COR34 = THETA (13)   ;Age-MSD
COR15 = THETA (14)   ;Dis-EXNB
COR25 = THETA (15)   ;slope-EXNB
COR35 = THETA (16)   ;Age-EXNB
COR45 = THETA (17)   ;MSD-EXNB


; Off-diagonal (Covariance)
COV12  = SQRT(OM11*OM22)*COR12
COV13  = SQRT(OM11*OM33)*COR13
COV14  = SQRT(OM11*OM44)*COR14
COV15  = SQRT(OM11*OM55)*COR15
COV23  = SQRT(OM22*OM33)*COR23
COV24  = SQRT(OM22*OM44)*COR24
COV25  = SQRT(OM22*OM55)*COR25
COV34  = SQRT(OM33*OM44)*COR34
COV35  = SQRT(OM33*OM55)*COR35
COV45  = SQRT(OM44*OM55)*COR45


; "Interaction terms"
C12   = COV12/SQRT(OM11)
C13   = COV13/SQRT(OM11)
C14   = COV14/SQRT(OM11)
C15   = COV15/SQRT(OM11)
C23   =(COV23-(C12*C13))/SQRT(OM22-(C12*C12))
C24   =(COV24-(C12*C14))/SQRT(OM22-(C12*C12))
C25   =(COV25-(C12*C15))/SQRT(OM22-(C12*C12))
C34   =(COV34-(C13*C14+C23*C24))/SQRT(OM33-(C13*C13+C23*C23))
C35   =(COV35-(C13*C15+C23*C25))/SQRT(OM33-(C13*C13+C23*C23))
C45   =(COV45-(C14*C15+C25*C24+C35*C34))/SQRT(OM44-(C14*C14+C24*C24+C34*C34))


P1     = SQRT(OM11)*ETA(1)                                                                                       ; Random effect Dis
P2     = C12*ETA(1) + SQRT(OM22-(C12*C12))*ETA(2)                                                                ; Random effect Slope
P3     = C13*ETA(1) + C23*ETA(2)+ SQRT(OM33-(C13*C13+C23*C23))*ETA(3)                                            ; Random effect Age
P4     = C14*ETA(1) + C24*ETA(2)+ C34*ETA(3)+ SQRT(OM44-(C14*C14+C24*C24+C34*C34))*ETA(4)                        ; Random effect MSD
P5     = C15*ETA(1) + C25*ETA(2)+ C35*ETA(3)+ C45*ETA(4)+ SQRT(OM55-(C15*C15+C25*C25+C35*C35+C45*C45))*ETA(5)    ; Random effect EXNB

; Drug effect
; Exposure dependent drug effect parameters
EFSM  = THETA(18)*EXP(ETA(6))
EC5S=THETA(19)

; Exposure independent drug effect
EFPM=THETA(20)

; Surrogate exposure calculations
CRL=CRCL
IF(CRCL.GT.150) CRL=150
EXPS=CD*104.5/CRL

EFSS=0
IF(TRT.GE.1.AND.DSS.GT.0) EFSS=EFSM*EXPS/(EXPS+EC5S)

EFPP=0
IF(TRT.GE.1.AND.DSS.GT.0) EFPP=EFPM

PD = P1 + ((THETA(1)+ P2)*(TIME/365)**(THETA(2)))*(1-EFPP)-EFSS            ;Disease parameter

;ICC parameters for EDSS subscores
IF(FLAGFS.EQ.1) THEN                   ;values item 1:PYRAMIDAL 0-5
BGE1   =THETA(24)
BGE2   =THETA(25)
BGE3   =THETA(26)
BGE4   =THETA(27)
BGE5   =THETA(28)
A_1  =THETA(29)
ENDIF
IF(FLAGFS.EQ.2) THEN                   ;values item 2:CEREBELLAR 0-5
BGE1   =THETA(30)
BGE2   =THETA(31)
BGE3   =THETA(32)
BGE4   =THETA(33)
BGE5   =THETA(34)
A_1  =THETA(35)
ENDIF

IF(FLAGFS.EQ.3) THEN                   ;values item 3:BRAINSTEM  0-4
BGE1   =THETA(36)
BGE2   =THETA(37)
BGE3   =THETA(38)
BGE4   =THETA(39)
A_1  =THETA(40)
ENDIF

IF(FLAGFS.EQ.4) THEN                   ;values item 4:SENSORY 0-6
BGE1   =THETA(41)
BGE2   =THETA(42)
BGE3   =THETA(43)
BGE4   =THETA(44)
BGE5   =THETA(45)
BGE6   =THETA(46)
A_1  =THETA(47)
ENDIF

IF(FLAGFS.EQ.5) THEN                   ;values item 5:BOWEL&BLADDER  0-5
BGE1   =THETA(48)
BGE2   =THETA(49)
BGE3   =THETA(50)
BGE4   =THETA(51)
BGE5   =THETA(52)
A_1  =THETA(53)
ENDIF

IF(FLAGFS.EQ.6) THEN                   ;values item 6: VISUAL  0-6
BGE1   =THETA(54)
BGE2   =THETA(55)
BGE3   =THETA(56)
BGE4   =THETA(57)
BGE5   =THETA(58)
BGE6   =THETA(59)
A_1  =THETA(60)
ENDIF

IF(FLAGFS.EQ.7) THEN                   ;values item 7:MENTAL  0-4
BGE1   =THETA(61)
BGE2   =THETA(62)
BGE3   =THETA(63)
BGE4   =THETA(64)
A_1  =THETA(65)
ENDIF

IF(FLAGFS.EQ.8) THEN                   ;values item 8:AMBAID 0-9
BGE1   =THETA(66)
BGE2   =THETA(67)
BGE3   =THETA(68)
BGE4   =THETA(69)
BGE5   =THETA(70)
BGE6   =THETA(71)
BGE7   =THETA(72)
BGE8   =THETA(73)
BGE9   =THETA(74)
A_1  =THETA(75)
ENDIF

Y=EPS(1)
IF(RTYPE.EQ.0) F_FLAG=2 ; -2 log likelihood

IF(FLAGFS.EQ.1) THEN        ;0-5
LGE1 =A_1*(PD-BGE1)                         ;Logits
LGE2 =A_1*(PD-BGE1-BGE2)
LGE3 =A_1*(PD-BGE1-BGE2-BGE3)
LGE4 =A_1*(PD-BGE1-BGE2-BGE3-BGE4)
LGE5 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5)

PGE1   =EXP(LGE1)/(1+EXP(LGE1))             ;probabilities
PGE2   =EXP(LGE2)/(1+EXP(LGE2))
PGE3   =EXP(LGE3)/(1+EXP(LGE3))
PGE4   =EXP(LGE4)/(1+EXP(LGE4))
PGE5   =EXP(LGE5)/(1+EXP(LGE5))

P01 = 1-PGE1
P11 = PGE1-PGE2
P21 = PGE2-PGE3
P31 = PGE3-PGE4
P41 = PGE4-PGE5
P51 = PGE5
ENDIF

P=1E-16
IF(FLAGFS.EQ.1.AND.DV.EQ.0) P=P01
IF(FLAGFS.EQ.1.AND.DV.EQ.1) P=P11
IF(FLAGFS.EQ.1.AND.DV.EQ.2) P=P21
IF(FLAGFS.EQ.1.AND.DV.EQ.3) P=P31
IF(FLAGFS.EQ.1.AND.DV.EQ.4) P=P41
IF(FLAGFS.EQ.1.AND.DV.GE.5) P=P51


IF(FLAGFS.EQ.2) THEN        ;0-5
LGE1 =A_1*(PD-BGE1)                         ;Logits
LGE2 =A_1*(PD-BGE1-BGE2)
LGE3 =A_1*(PD-BGE1-BGE2-BGE3)
LGE4 =A_1*(PD-BGE1-BGE2-BGE3-BGE4)
LGE5 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5)

PGE1   =EXP(LGE1)/(1+EXP(LGE1))             ;probabilities
PGE2   =EXP(LGE2)/(1+EXP(LGE2))
PGE3   =EXP(LGE3)/(1+EXP(LGE3))
PGE4   =EXP(LGE4)/(1+EXP(LGE4))
PGE5   =EXP(LGE5)/(1+EXP(LGE5))

P02 = 1-PGE1
P12 = PGE1-PGE2
P22 = PGE2-PGE3
P32 = PGE3-PGE4
P42 = PGE4-PGE5
P52 = PGE5
ENDIF

IF(FLAGFS.EQ.2.AND.DV.EQ.0) P=P02
IF(FLAGFS.EQ.2.AND.DV.EQ.1) P=P12
IF(FLAGFS.EQ.2.AND.DV.EQ.2) P=P22
IF(FLAGFS.EQ.2.AND.DV.EQ.3) P=P32
IF(FLAGFS.EQ.2.AND.DV.EQ.4) P=P42
IF(FLAGFS.EQ.2.AND.DV.GE.5) P=P52


IF(FLAGFS.EQ.3) THEN        ;0-4
LGE1 =A_1*(PD-BGE1)                         ;Logits
LGE2 =A_1*(PD-BGE1-BGE2)
LGE3 =A_1*(PD-BGE1-BGE2-BGE3)
LGE4 =A_1*(PD-BGE1-BGE2-BGE3-BGE4)

PGE1   =EXP(LGE1)/(1+EXP(LGE1))             ;probabilities
PGE2   =EXP(LGE2)/(1+EXP(LGE2))
PGE3   =EXP(LGE3)/(1+EXP(LGE3))
PGE4   =EXP(LGE4)/(1+EXP(LGE4))

P03 = 1-PGE1
P13 = PGE1-PGE2
P23 = PGE2-PGE3
P33 = PGE3-PGE4
P43 = PGE4
ENDIF

IF(FLAGFS.EQ.3.AND.DV.EQ.0) P=P03
IF(FLAGFS.EQ.3.AND.DV.EQ.1) P=P13
IF(FLAGFS.EQ.3.AND.DV.EQ.2) P=P23
IF(FLAGFS.EQ.3.AND.DV.EQ.3) P=P33
IF(FLAGFS.EQ.3.AND.DV.GE.4) P=P43


IF(FLAGFS.EQ.4) THEN        ;0-6
LGE1 =A_1*(PD-BGE1)                         ;Logits
LGE2 =A_1*(PD-BGE1-BGE2)
LGE3 =A_1*(PD-BGE1-BGE2-BGE3)
LGE4 =A_1*(PD-BGE1-BGE2-BGE3-BGE4)
LGE5 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5)
LGE6 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5-BGE6)

PGE1   =EXP(LGE1)/(1+EXP(LGE1))             ;probabilities
PGE2   =EXP(LGE2)/(1+EXP(LGE2))
PGE3   =EXP(LGE3)/(1+EXP(LGE3))
PGE4   =EXP(LGE4)/(1+EXP(LGE4))
PGE5   =EXP(LGE5)/(1+EXP(LGE5))
PGE6   =EXP(LGE6)/(1+EXP(LGE6))

P04 = 1-PGE1
P14 = PGE1-PGE2
P24 = PGE2-PGE3
P34 = PGE3-PGE4
P44 = PGE4-PGE5
P54 = PGE5-PGE6
P64 = PGE6
ENDIF

IF(FLAGFS.EQ.4.AND.DV.EQ.0) P=P04
IF(FLAGFS.EQ.4.AND.DV.EQ.1) P=P14
IF(FLAGFS.EQ.4.AND.DV.EQ.2) P=P24
IF(FLAGFS.EQ.4.AND.DV.EQ.3) P=P34
IF(FLAGFS.EQ.4.AND.DV.EQ.4) P=P44
IF(FLAGFS.EQ.4.AND.DV.EQ.5) P=P54
IF(FLAGFS.EQ.4.AND.DV.EQ.6) P=P64


IF(FLAGFS.EQ.5) THEN        ;0-5
LGE1 =A_1*(PD-BGE1)                         ;Logits
LGE2 =A_1*(PD-BGE1-BGE2)
LGE3 =A_1*(PD-BGE1-BGE2-BGE3)
LGE4 =A_1*(PD-BGE1-BGE2-BGE3-BGE4)
LGE5 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5)

PGE1   =EXP(LGE1)/(1+EXP(LGE1))             ;probabilities
PGE2   =EXP(LGE2)/(1+EXP(LGE2))
PGE3   =EXP(LGE3)/(1+EXP(LGE3))
PGE4   =EXP(LGE4)/(1+EXP(LGE4))
PGE5   =EXP(LGE5)/(1+EXP(LGE5))

P05 = 1-PGE1
P15 = PGE1-PGE2
P25 = PGE2-PGE3
P35 = PGE3-PGE4
P45 = PGE4-PGE5
P55 = PGE5
ENDIF

IF(FLAGFS.EQ.5.AND.DV.EQ.0) P=P05
IF(FLAGFS.EQ.5.AND.DV.EQ.1) P=P15
IF(FLAGFS.EQ.5.AND.DV.EQ.2) P=P25
IF(FLAGFS.EQ.5.AND.DV.EQ.3) P=P35
IF(FLAGFS.EQ.5.AND.DV.EQ.4) P=P45
IF(FLAGFS.EQ.5.AND.DV.GE.5) P=P55


IF(FLAGFS.EQ.6) THEN        ;0-6
LGE1 =A_1*(PD-BGE1)                         ;Logits
LGE2 =A_1*(PD-BGE1-BGE2)
LGE3 =A_1*(PD-BGE1-BGE2-BGE3)
LGE4 =A_1*(PD-BGE1-BGE2-BGE3-BGE4)
LGE5 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5)
LGE6 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5-BGE6)

PGE1   =EXP(LGE1)/(1+EXP(LGE1))             ;probabilities
PGE2   =EXP(LGE2)/(1+EXP(LGE2))
PGE3   =EXP(LGE3)/(1+EXP(LGE3))
PGE4   =EXP(LGE4)/(1+EXP(LGE4))
PGE5   =EXP(LGE5)/(1+EXP(LGE5))
PGE6   =EXP(LGE6)/(1+EXP(LGE6))

P06 = 1-PGE1
P16 = PGE1-PGE2
P26 = PGE2-PGE3
P36 = PGE3-PGE4
P46 = PGE4-PGE5
P56 = PGE5-PGE6
P66 = PGE6
ENDIF

IF(FLAGFS.EQ.6.AND.DV.EQ.0) P=P06
IF(FLAGFS.EQ.6.AND.DV.EQ.1) P=P16
IF(FLAGFS.EQ.6.AND.DV.EQ.2) P=P26
IF(FLAGFS.EQ.6.AND.DV.EQ.3) P=P36
IF(FLAGFS.EQ.6.AND.DV.EQ.4) P=P46
IF(FLAGFS.EQ.6.AND.DV.EQ.5) P=P56
IF(FLAGFS.EQ.6.AND.DV.EQ.6) P=P66


IF(FLAGFS.EQ.7) THEN        ;0-4
LGE1 =A_1*(PD-BGE1)                         ;Logits
LGE2 =A_1*(PD-BGE1-BGE2)
LGE3 =A_1*(PD-BGE1-BGE2-BGE3)
LGE4 =A_1*(PD-BGE1-BGE2-BGE3-BGE4)

PGE1   =EXP(LGE1)/(1+EXP(LGE1))             ;probabilities
PGE2   =EXP(LGE2)/(1+EXP(LGE2))
PGE3   =EXP(LGE3)/(1+EXP(LGE3))
PGE4   =EXP(LGE4)/(1+EXP(LGE4))

P07 = 1-PGE1
P17 = PGE1-PGE2
P27 = PGE2-PGE3
P37 = PGE3-PGE4
P47 = PGE4
ENDIF

IF(FLAGFS.EQ.7.AND.DV.EQ.0) P=P07
IF(FLAGFS.EQ.7.AND.DV.EQ.1) P=P17
IF(FLAGFS.EQ.7.AND.DV.EQ.2) P=P27
IF(FLAGFS.EQ.7.AND.DV.EQ.3) P=P37
IF(FLAGFS.EQ.7.AND.DV.GE.4) P=P47


IF(FLAGFS.EQ.8) THEN        ;0-9
LGE1 =A_1*(PD-BGE1)                         ;Logits
LGE2 =A_1*(PD-BGE1-BGE2)
LGE3 =A_1*(PD-BGE1-BGE2-BGE3)
LGE4 =A_1*(PD-BGE1-BGE2-BGE3-BGE4)
LGE5 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5)
LGE6 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5-BGE6)
LGE7 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5-BGE6-BGE7)
LGE8 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5-BGE6-BGE7-BGE8)
LGE9 =A_1*(PD-BGE1-BGE2-BGE3-BGE4-BGE5-BGE6-BGE7-BGE8-BGE9)

PGE1   =EXP(LGE1)/(1+EXP(LGE1))             ;probabilities
PGE2   =EXP(LGE2)/(1+EXP(LGE2))
PGE3   =EXP(LGE3)/(1+EXP(LGE3))
PGE4   =EXP(LGE4)/(1+EXP(LGE4))
PGE5   =EXP(LGE5)/(1+EXP(LGE5))
PGE6   =EXP(LGE6)/(1+EXP(LGE6))
PGE7   =EXP(LGE7)/(1+EXP(LGE7))
PGE8   =EXP(LGE8)/(1+EXP(LGE8))
PGE9   =EXP(LGE9)/(1+EXP(LGE9))

P08 = 1-PGE1
P18 = PGE1-PGE2
P28 = PGE2-PGE3
P38 = PGE3-PGE4
P48 = PGE4-PGE5
P58 = PGE5-PGE6
P68 = PGE6-PGE7
P78 = PGE7-PGE8
P88 = PGE8-PGE9
P98 = PGE9
ENDIF

IF(FLAGFS.EQ.8.AND.DV.EQ.0) P=P08
IF(FLAGFS.EQ.8.AND.DV.EQ.1) P=P18
IF(FLAGFS.EQ.8.AND.DV.EQ.2) P=P28
IF(FLAGFS.EQ.8.AND.DV.EQ.3) P=P38
IF(FLAGFS.EQ.8.AND.DV.EQ.4) P=P48
IF(FLAGFS.EQ.8.AND.DV.EQ.5) P=P58
IF(FLAGFS.EQ.8.AND.DV.EQ.6) P=P68
IF(FLAGFS.EQ.8.AND.DV.EQ.7) P=P78
IF(FLAGFS.EQ.8.AND.DV.EQ.8) P=P88
IF(FLAGFS.EQ.8.AND.DV.GE.9) P=P98

IF(P.LT.1E-16) P = 1E-16      ;protection against low values
Y=-2*LOG(P)

;FREM
IF(RTYPE.EQ.1) Y=THETA(21)+ P3 +EPS(1)      ; Age
IF(RTYPE.EQ.2) Y=THETA(22)+ P4 +EPS(1)      ; MSD
IF(RTYPE.EQ.3) Y=THETA(23)+ P5 +EPS(1)      ; EXNB

;Simulations
;IF (ICALL.EQ.4) THEN
;    CALL RANDOM (2,R)
;   RAND=R
;   SIM=1
;    IF(R.LE.PGE2) SIM=2
;    IF(R.LE.PGE3) SIM=3
;    IF(R.LE.PGE4) SIM=4
;    IF(R.LE.PGE5) SIM=5
;    IF(R.LE.PGE6) SIM=6
;ENDIF


$THETA
(0, 0.093)       ;1 disease prog slope
(0, 0.710)       ;2 disease prog power
(1) FIX          ;3 var_P1 Disability
(0, 0.2)         ;4 var_P2 Slope
(0,99.329)       ;5 Var_P3 AGE
(0,54.124)       ;6 Var_P4 MSD
(0,0.3727)       ;7 Var_P5 EXNB
(-1, 0.11,1)     ;8 Cor Dis-slope
(-1, 0.2645,1)   ;9 Cor Dis-age
(-1, 0.1175,1)   ;10 Cor slope-Age
(-1, 0.2727,1)   ;11 Cor Dis-MSD
(-1, 0.0896,1)   ;12 Cor slope-MSD
(-1, 0.4582,1)   ;13 Cor age-MSD
(-1, 0.0391,1)   ;14 COR dis-EXNB
(-1, 0.0695,1)   ;15 COR slope-EXNB
(-1, -0.0914,1)  ;16 COR age-EXNB
(-1, -0.115,1)   ;17 Cor MSD-EXNB
(0, 0.17)        ;18 Emax for symptomatic effect
(0, 408.29)      ;19 EC50 for symptomatic effect
(0, 0.209,1)     ;20 protective effecT
(0, 38.569)      ;21 Age Mean
(0, 8.6697)      ;22 MSD Mean
(0, 1.3540)      ;23 EXNB Mean
(-1.55)          ;24 item 1b1
(0,1.248)        ;25 item 1b2
(0,0.818)        ;26 item 1b3
(0,1.428)        ;27 item 1b4
(0,1.313)        ;28 item 1b5
(0,3.172)        ;29 item 1slope
(-0.913)         ;30 item 2b1
(0,0.963)        ;31 item 2b2
(0,1.023)        ;32 item 2b3
(0,1.648)        ;33 item 2b4
(0,1.074)        ;34 item 2b5
(0,2.873)        ;35 item 2slope
(-0.112)         ;36 item 3b1
(0,1.711)        ;37 item 3b2
(0,2.000)        ;38 item 3b3
(0,2.780)        ;39 item 3b4
(0,1.038)        ;40 item 3slope
(-0.783)         ;41 item 4b1
(0,1.186)        ;42 item 4b2
(0,1.946)        ;43 item 4b3
(0,2.345)        ;44 item 4b4
(0,2.937)        ;45 item 4b5
(0,2.353)        ;46 item 4b6
(0,0.993)        ;47 item 4slope
(-0.147)         ;48 item 5b1
(0,1.887)        ;49 item 5b2
(0,1.531)        ;50 item 5b3
(0,2.996)        ;51 item 5b4
(0,0.968)        ;52 item 5b5
(0,1.256)        ;53 item 5slope
(-0.037)         ;54 item 6b1
(0,3.751)        ;55 item 6b2
(0,2.660)        ;56 item 6b3
(0,1.6877)       ;57 item 6b4
(0,1.705)        ;58 item 6b5
(0,1.637)        ;59 item 6b6
(0,0.440)        ;60 item 6slope
(0.402)          ;61 item 7b1
(0,1.111)        ;62 item 7b2
(0,4.161)        ;63 item 7b3
(0,2.890)        ;64 item 7b4
(0,0.912)        ;65 item 7slope
(1.133)          ;66 item 8b1
(0,0.244)        ;67 item 8b2
(0,0.225)        ;68 item 8b3
(0,0.423)        ;69 item 8b4
(0,0.448)        ;70 item 8b5
(0,0.484)        ;71 item 8b6
(0,0.139)        ;72 item 8b7
(0,0.261)        ;73 item 8b8
(0,0.351)        ;74 item 8b9
(0,3.64)         ;75 item 8slope

$OMEGA 1 FIX     ;1 Disability
$OMEGA 1 FIX     ;2 Slope
$OMEGA 1 FIX     ;3 Age
$OMEGA 1 FIX     ;4 MSD
$OMEGA 1 FIX     ;5 EXNB
$OMEGA 2.06      ;4 Emax

$SIGMA 0.00001 FIX


$ESTIM MAXEVAL=0 METHOD=COND LAPLACE PRINT=10 NOABORT SLOW MSFO=msfb1002h ;
;$COV MATRIX=R
;$SIMULATION (1707)(071981 UNIFORM) ONLYSIM NSUB=1
;$TABLE ID TIME VISIT EDSS TRT PD FLAGFS DV NOPRINT NOAPPEND ONEHEADER FILE=sdtab1002h