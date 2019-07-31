$PROB   Prediction of daily ANC

$INPUT  NREF        ; Record identification number
        ID          ; ID number
        DV          ; Simulated DV (log scale)
        TIME        ; hours
        CP          ; Predicted docetaxel cocnetration
        EV21=EVID   ; EVID to use for making prediction based on monitoring up to 21 days
        EV20=DROP   ; EVID to use for making prediction based on monitoring up to 20 days
        EV19=DROP   ; EVID to use for making prediction based on monitoring up to 19 days
        EV18=DROP   ; EVID to use for making prediction based on monitoring up to 18 days
        EV17=DROP   ; EVID to use for making prediction based on monitoring up to 17 days
        EV16=DROP   ; EVID to use for making prediction based on monitoring up to 16 days
        EV15=DROP   ; EVID to use for making prediction based on monitoring up to 15 days
        EV14=DROP   ; EVID to use for making prediction based on monitoring up to 14 days
        EV13=DROP   ; EVID to use for making prediction based on monitoring up to 13 days
        EV12=DROP   ; EVID to use for making prediction based on monitoring up to 12 days
        EV11=DROP   ; EVID to use for making prediction based on monitoring up to 11 days
        EV10=DROP   ; EVID to use for making prediction based on monitoring up to 10 days
        EV9=DROP    ; EVID to use for making prediction based on monitoring up to 9 days
        EV8=DROP    ; EVID to use for making prediction based on monitoring up to 8 days
        EV7=DROP    ; EVID to use for making prediction based on monitoring up to 7 days
        EV6=DROP    ; EVID to use for making prediction based on monitoring up to 6 days
        EV5=DROP    ; EVID to use for making prediction based on monitoring up to 5 days
        EV4=DROP    ; EVID to use for making prediction based on monitoring up to 4 days
        EV3=DROP    ; EVID to use for making prediction based on monitoring up to 3 days
        DAY         ; Dummy variable (1 if time for prediction, 0 if other)
        AMT         ; Dummy amount
        CMT
        SEX         ;
        PERF        ; Performance status (categorical)
        PC          ; Previous anticancer therapy (categorical)
        AAG         ; alpha1-acid glycoprotein (continuous)

$DATA Simulated_myelosuppression_dailyANC.csv
      IGN=@

$SUBS   ADVAN6 TOL=5
$MODEL  NCOMP=12
        COMP=(CIRC,DEFOBS) ; Circulating neutrophils
        COMP=(STEM)        ; Proliferating stem cells
        COMP=(TRANSIT1)    ; Cells in transit compartment 1
        COMP=(TRANSIT2)    ; Cells in transit compartment 2
        COMP=(TRANSIT3)    ; Cells in transit compartment 3
        COMP=(RTBA)        ; Return to baseline
        COMP=(OOG0)        ; Occurence of Grade 0
        COMP=(OOG4)        ; Occurence of Grade 4
        COMP=(EOG4)        ; End of Grade 4
        COMP=(OO01)        ; Occurence of ANC 0.1
        COMP=(ONADIR)      ; Occurence of NADIR
        COMP=(NADIR)       ; NADIR value

$PK
"FIRST
" COMMON /PRCOMG/ IDUM1,IDUM2,IMAX,IDUM4,IDUM5
" INTEGER IDUM1,IDUM2,IMAX,IDUM4,IDUM5
" IMAX=70000000

; =====================================================
; Covariate effects as according to Kloft et al., 2006
; =====================================================

IF(SEX.EQ.1) BASEX = 0
IF(SEX.EQ.2) BASEX = THETA(11)

IF(PERF.EQ.0.OR.PERF.EQ.-99) BAPERF = 0
IF(PERF.GE.1) BAPERF = THETA(8)

IF(PC.EQ.0) BAPC = 0
IF(PC.EQ.1) BAPC = THETA(7)

IF(AAG.LE.1.34) BAAAG = THETA(9)*(AAG-1.34)
IF(AAG.GT.1.34) BAAAG = THETA(10)*(AAG-1.34)

SLAAG = THETA(6)*(AAG-1.34)

BACOV=(1+BASEX)*(1+BAPERF)*(1+BAPC)
BACOV=BACOV*(1+BAAAG)

MTCOV=1

SLCOV=(1+SLAAG)

; =====================================================
; Parameters in the myelosuppression model
; =====================================================

BA = THETA(1)*EXP(ETA(1))*BACOV
MT = THETA(2)*EXP(ETA(2))*MTCOV
SL = THETA(3)/808*1000*EXP(ETA(3))*SLCOV
PO = THETA(4)

K  =  (4/MT)
F1 =  BA
F2 =  BA
F3 =  BA
F4 =  BA
F5 =  BA

$DES

; =====================================================
; Myelosuppression model
; =====================================================

; Linear drug effect
DRUG = SL*CP

D1       =  K*A(5) - K*A(1)
DADT(1)  =  D1                                     ; circulating neutrophils
DADT(2)  =  K*A(2)*(1-DRUG)*(BA/A(1))**PO - K*A(2) ; proliferating cells
DADT(3)  =  K*A(2) - K*A(3)                        ; transit 1
DADT(4)  =  K*A(3) - K*A(4)                        ; transit 2
DADT(5)  =  K*A(4) - K*A(5)                        ; transit 3

; =================================================================
; Compute markers of summary variables and other neutropenic events
; =================================================================

NEUL = A(1)
B    = K*A(5)-K*A(1)
                                                   ; "slope" used to extract model-predicted nadir
; Marks return to baseline
IF(NEUL.GT.BA)THEN
  DADT(6)  =  1
ELSE
  DADT(6)  =  0
ENDIF

; Marks occurence of Grade 0 (ANC=2)
IF(NEUL.GE.2.AND.B.GT.0)THEN
  DADT(7) = 1
ELSE
  DADT(7) = 0
ENDIF

; Marks occurence of Grade 4 (ANC=0.5)
IF(NEUL.LT.0.5)THEN
  DADT(8) = 1
ELSE
  DADT(8) = 0
ENDIF

; Marks end of Grade 4
IF(NEUL.GT.0.5.AND.B.GT.0)THEN
  DADT(9) = 1
ELSE
  DADT(9) = 0
ENDIF

; Marks occurence of ANC<=0.1
IF(NEUL.LE.0.1)THEN
  DADT(10) = 1
ELSE
  DADT(10) = 0
ENDIF

; Marks nadir
IF(B.GE.0.AND.TIME.GT.72)THEN
  DADT(11) = 1                                    ; "Time to nadir"
  DADT(12) = D1                                   ; "nadir value"
ELSE
  DADT(11) = 0
  DADT(12) = 0
ENDIF


$ERROR
  DEL    = 0
  IF(F.EQ.0)DEL = 1
  IPRED  =  LOG(F+DEL)
  W      =  THETA(5)
  IRES   =  DV -  IPRED
  IWRES  =  IRES/W
  Y      = LOG(F)+W*EPS(1)

AA1    = A(1)  ; circulating neutrophils
AA2    = A(2)  ; proliferating cells
AA3    = A(3)  ; transit 1
AA4    = A(4)  ; transit 2
AA5    = A(5)  ; transit 3
RTBA   = A(6)  ; Marker for return to baseline
OOG0   = A(7)  ; Marker for occurence of grade 0 neutropenia
OOG4   = A(8)  ; Marker for occurence of grade 4 neutropenia
EOG4   = A(9)  ; Marker for end of grade 4 neutropenia
OO01   = A(10) ; Marker for occurence of ANC<=0.1
ONADIR = A(11) ; Marker for nadir (time)
NADIR  = A(12) ; Marker for nadir (ANC)

; =====================================================
; Computing values of summary variables
; =====================================================

; NB: in the code below LE to 24.001 is used because how the data is set up
FLG6=0
IF(A(6).GT.0.AND.A(6).LE.24.001) FLG6=1
TTBA=0                                     ; Time to baseline
IF(FLG6.EQ.1) TTBA = TIME - A(6)

FLG7=0
IF(A(7).GT.0.AND.A(7).LE.24.001) FLG7=1
TTG0=0                                     ; Time to grade 0
IF(FLG7.EQ.1) TTG0 = TIME - A(7)

FLG8=0
IF(A(1).LT.0.5.AND.A(8).LE.24.001.AND.EVID.EQ.0) FLG8=1
TTG4=0                                     ; Time to grade 4
IF(FLG8.EQ.1) TTG4 = TIME - A(8)

FLG9=0
IF(A(9).GT.0.AND.A(9).LE.24.001) FLG9=1
DOG4=0                                     ; Duration of grade 4
IF(FLG9.EQ.1) DOG4 = OOG4

FLG10=0
IF(A(1).LT.0.1.AND.A(10).GT.0.AND.A(10).LE.24.001) FLG10=1
TT01=0                                      ; Time to ANC 0.1
IF(FLG10.EQ.1) TT01 = TIME - A(10)

FLG11=0
IF(A(11).GT.0.AND.A(11).LE.24.001) FLG11=1
TNADIR=0                                    ; Time to nadir
VNADIR=0                                    ; nadir value
IF(FLG11.EQ.1) TNADIR = TIME - A(11)
IF(FLG11.EQ.1) VNADIR = A(1) - A(12)

; =====================================================
; Parameter estimates as according to Kloft et al., 2006
; =====================================================
$THETA  (3,5.21965)              ; 1 BA
$THETA  (60,84.1791)             ; 2 MT
$THETA  (0,15.5711)              ; 3 SL
$THETA  (0,0.144543)             ; 4 PO
$THETA  (0,0.424093)             ; 5 res err - prop on norm sc
$THETA  (-0.436,-0.350693,0.952) ; 6 SLAAG
$THETA  (-0.999,-0.146837,99)    ; 7 BAPC
$THETA  (-0.999,0.130406,99)     ; 8 BAPERF
$THETA  (-0.436,0.174677,0.952)  ; 9  BAAAG LE medAAG
$THETA  (-0.436,0.494618,0.952)  ; 10 BAAAG GT medAAG
$THETA  (-0.999,-0.121451,99)    ; 11 BASEX
$OMEGA  0.0639703                ; 1 IIV BASE
$OMEGA  0.0191785                ; 2 IIV MTT
$OMEGA  0.128412                 ; 3 IIV SLOPE
$SIGMA  1  FIX

$EST NOABORT METH=1 INTER SIGDIGITS=2 MAXEVALS=0
 PRINT=5 POSTHOC

$TABLE ID TIME DV EVID B AA1 RTBA FLG6 TTBA OOG0 FLG7
 TTG0 OOG4 FLG8 TTG4 EOG4 FLG9 DOG4 OO01 FLG10 TT01
 ONADIR NADIR FLG11 TNADIR VNADIR NOPRINT ONEHEADER FILE=Predicted_myelosuppression_dailyANC

