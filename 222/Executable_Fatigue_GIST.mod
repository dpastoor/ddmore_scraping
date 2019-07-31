$PROBLEM Published model CPT_PSP

$INPUT ID  TIME DV MDV DOSE CL BAS3 MRT3 EC53

$DATA     Simulated_Fatigue_GIST.txt  IGNORE=#

$SUBROUTINE ADVAN6 TOL=6

$MODEL       NCOMPS=1
             COMP=COMP1 ;VEGF3

$PK

  A_0(1)   =         BAS3

  EMAX     =         1
  HILL     =         3.31
  HILL2    =         1.54

  KOUT3    =         1/MRT3
  KIN3     =         BAS3*KOUT3

  AUC      =         DOSE/CL


$DES

;-------------------PD-----------------------

   EFF3     =      EMAX*AUC/(EC53+AUC)
   DADT(1)  =      KIN3*(1-EFF3)-KOUT3*A(1)     ;VEGFR3

$ERROR

;-------------------BIOMARKER----------------------
    VE3   =  A(1)                     ; sVEGFR-3
    BM    = ((VE3-BAS3)/BAS3)         ; Relative change in sVEGFR-3

;-------BIOMARKER
   VE30 = THETA(13)*BM
   VE31 = THETA(14)*BM
   VE32 = THETA(15)*BM
   VE33 = THETA(16)*BM

;--------INIT----------------------------------------------
 IF(NEWIND.NE.2) TMP=DV    ; Observation at time=0
      PDV=TMP              ; Previous score
;---------------------------------------------------------------

; PX0        Transition from 0 to 1, 0 to 2 or 0 to 3
   B10 = THETA(1) + VE30 + ETA(1)       ; Predictor of PX0 Transitions
   B20 = THETA(2)
   B30 = THETA(3)

   LGE10 = B10
   LGE20 = B10 + B20
   LGE30 = B10 + B20 + B30

;Logits
   PGE10  = EXP(LGE10)/(1+EXP(LGE10))
   PGE20  = EXP(LGE20)/(1+EXP(LGE20))
   PGE30  = EXP(LGE30)/(1+EXP(LGE30))


; PX1        Transition from 1 to 1, 1 to 2 or 1 to 3
   B11 = THETA(4) + VE31 + ETA(2)       ; Predictor of PX1 Transitions
   B21 = THETA(5)
   B31 = THETA(6)

   LGE11 = B11
   LGE21 = B11 + B21
   LGE31 = B11 + B21 + B31

;Logits
   PGE11  = EXP(LGE11)/(1+EXP(LGE11))
   PGE21  = EXP(LGE21)/(1+EXP(LGE21))
   PGE31  = EXP(LGE31)/(1+EXP(LGE31))

; PX2        Transition from 2 to 1, 2 to 2 or 2 to 3
   B12 = THETA(7) + VE32 + ETA(3)       ; Predictor of PX2 Transitions
   B22 = THETA(8)
   B32 = THETA(9)

   LGE12 = B12
   LGE22 = B12 + B22
   LGE32 = B12 + B22 + B32

;Logits
   PGE12  = EXP(LGE12)/(1+EXP(LGE12))
   PGE22  = EXP(LGE22)/(1+EXP(LGE22))
   PGE32  = EXP(LGE32)/(1+EXP(LGE32))


; PX3        Transition from 3 to 1, 3 to 2 or 3 to 3
   B13 = THETA(10) + VE33 + ETA(4)      ; Predictor of PX3 Transitions
   B23 = THETA(11)
   B33 = THETA(12)

   LGE13 = B13
   LGE23 = B13 + B23
   LGE33 = B13 + B23 + B33

;Logits
   PGE13  = EXP(LGE13)/(1+EXP(LGE13))
   PGE23  = EXP(LGE23)/(1+EXP(LGE23))
   PGE33  = EXP(LGE33)/(1+EXP(LGE33))
   
;=====================================
   P00 = (1-PGE10)
   P10 = (PGE10-PGE20)
   P20 = (PGE20-PGE30)
   P30 = PGE30
   
   P01 = (1-PGE11)
   P11 = (PGE11-PGE21)
   P21 = (PGE21-PGE31)
   P31 = PGE31

   P02 = (1-PGE12)
   P12 = (PGE12-PGE22)
   P22 = (PGE22-PGE32)
   P32 =  PGE32
   
   P03 = (1-PGE13)
   P13 = (PGE13-PGE23)
   P23 = (PGE23-PGE33)
   P33 = PGE33

;----PX0
IF(PDV.EQ.0.AND.DV.EQ.1) Y=P10
IF(PDV.EQ.0.AND.DV.EQ.2) Y=P20
IF(PDV.EQ.0.AND.DV.GT.2) Y=P30
IF(PDV.EQ.0.AND.DV.EQ.0) Y=P00

;----PX1
IF(PDV.EQ.1.AND.DV.EQ.0) Y=P01
IF(PDV.EQ.1.AND.DV.EQ.2) Y=P21
IF(PDV.EQ.1.AND.DV.GT.2) Y=P31
IF(PDV.EQ.1.AND.DV.EQ.1) Y=P11

;----PX2
IF(PDV.EQ.2.AND.DV.EQ.0) Y=P02
IF(PDV.EQ.2.AND.DV.EQ.1) Y=P12
IF(PDV.EQ.2.AND.DV.GT.2) Y=P32
IF(PDV.EQ.2.AND.DV.EQ.2) Y=P22

;----PX3
IF(PDV.GT.2.AND.DV.EQ.0) Y=P03
IF(PDV.GT.2.AND.DV.EQ.1) Y=P13
IF(PDV.GT.2.AND.DV.EQ.2) Y=P23
IF(PDV.GT.2.AND.DV.GT.2) Y=P33

;--------
        TMP=DV         ; To remember previous DV

$THETA
 -5.61000              ; 1  PX0 B1
 (-1000000,-1.14000,0) ; 2  PX0 B2
 (-1000000,-1.60000,0) ; 3  PX0 B3
 -2.78000              ; 4  PX1 B1
 (-1000000,-1.76000,0) ; 5  PX1 B2
 (-1000000,-1.77000,0) ; 6  PX1 B3
 -3.30000              ; 7  PX2 B1
 (-1000000,-0.96100,0) ; 8  PX2 B2
 (-1000000,-1.45000,0) ; 9  PX2 B3
 -3.38000              ; 10 PX3 B1
 (-1000000,-0.71600,0) ; 11 PX3 B2
 (-1000000,-0.09740,0) ; 12 PX3 B3


 .01                   ; 13 B10
 .01                   ; 14 B11
 .01                   ; 15 B12
 .01                   ; 16 B13


$OMEGA .01   ; 1 PX0
$OMEGA .01   ; 2 PX1
$OMEGA .01   ; 3 PX2
$OMEGA .01   ; 4 PX3


$ESTIMATION MAXEVAL=0 METHOD=1 LAPLACE LIKE PRINT=1 NOABORT
