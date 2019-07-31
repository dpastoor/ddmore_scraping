$SIZES      NO=800 LIM6=800     ; needed for big datasets
$PROBLEM Biomarker models for VEGF, VEGFR2, VEGFR3 and SKIT

; Reference: Hansson E.K. 2011. Pharmacometric models for Biomarkers, Side Effects and Efficacy in Anticancer Drug Therapy. Acta Universitatis Upsaliensis.

$INPUT
ID           ; Patient identification
CYCL         ; Cycle number
TIME         ; Time in hours
WEEK         ; Time in weeks
FLAG         ; 1. Dose ; 4. Tumor size (SLD) (ignored) ; 5.VEGF ; 6.VEGFR2 ; 7.VEGFR3 ; 8. SKIT
DVX          ; DV on a linear scale
DV           ; log-transformed DV
DOS          ; Dose in mg
PLA          ; Placebo: 1. untreated, 0. treated
CL           ; posthoc total plasma clearance from a previously developed 2-compartment model
EVID         ; 0.PD assessment ; 2. other event


$DATA   data_ddmore_Biomarker_GIST_simulated.csv
        IGN=# ACCEPT(FLAG.EQ.1,FLAG.EQ.5,FLAG.EQ.6,FLAG.EQ.7,FLAG.EQ.8)

$SUBROUTINE ADVAN6 TOL=5
$ABBR DERIV2=NOCOMMON

$MODEL
 NCOMPS=4
 COMP=COMP1  ;VEGF
 COMP=COMP2  ;sVEGFR-2
 COMP=COMP3  ;sVEGFR-3
 COMP=COMP4  ;SKIT

$PK
; Verbatim code : changes the iteration maximum (IMAX) (default value 100000)
"FIRST
" COMMON /PRCOMG/ IDUM1,IDUM2,IMAX,IDUM4,IDUM5
" INTEGER IDUM1,IDUM2,IMAX,IDUM4,IDUM5
" IMAX=70000000

;------VEGF------------------------------------

  BM0    =   THETA(1)*EXP(ETA(1))          ; Baseline VEGF
  MRT    =   THETA(2)*EXP(ETA(5))
  IMAX1  =   THETA(3)
  IC50   =   THETA(4)*EXP(ETA(7))          ; Common IC50 for the four biomarkers which were found to be highly correlated
  HILL   =   THETA(5)                      ; Hill coefficient for VEGF sigmoid Imax model
  TVSLO  =   THETA(6)/1000
  DPSLO  =   TVSLO *EXP(ETA(11))           ; Slope for the disease progression model
;------VEGFR2----------------------------------
  BM02   =   THETA(7)*EXP(ETA(2))          ; Baseline VEGFR2
  MRT2   =   THETA(8)*EXP(ETA(5))
  IMAX2  =   THETA(3)
  IC502  =   THETA(4)*EXP(ETA(8))
  HILL2  =   THETA(9)                      ; Hill coefficient for VEGFR2 sigmoid IMAX model

;------VEGFR3----------------------------------
  BM03   =   THETA(10)*EXP(ETA(3))         ; Baseline VEGFR3
  MRT3   =   THETA(11)*EXP(ETA(5))
  IMAX3  =   THETA(3)
  IC503  =   THETA(4)*EXP(ETA(9))

;------SKIT------------------------------------

  BM0S   =   THETA(12)*EXP(ETA(4))         ; Baseline SKIT
  MRTS   =   THETA(13)*EXP(ETA(6))
  IMAXS  =   THETA(3)
  IC50S  =   THETA(4)*EXP(ETA(10))
  TVSLOS =   THETA(6)/1000
  DPSLOS =   TVSLOS *EXP(ETA(12))          ; Slope for the disease progression model

;------Compartment initialization--------------
  A_0(1) =   BM0                           ; VEGF
  A_0(2) =   BM02                          ; sVEGFR-2
  A_0(3) =   BM03                          ; sVEGFR-3
  A_0(4) =   BM0S                          ; sKIT



  KOUT   =   1/MRT
  KOUT2  =   1/MRT2
  KOUT3  =   1/MRT3
  KOUTS  =   1/MRTS

  KIN2   =   BM02*KOUT2
  KIN3   =   BM03*KOUT3

  AUC    =   DOS/CL
  DP1     =  BM0*(1+DPSLO*TIME)            ; Linear disease progression model for VEGF
  DPS    =   BM0S*(1+DPSLOS*TIME)          ; Linear disease progression model for SKIT


$DES

 ;-----PD: Sunitinib effect--------------------

 KIN     =   DP1*KOUT
 KINS    =   DPS*KOUTS

 EFF     =   IMAX1*AUC**HILL /(IC50**HILL+AUC**HILL)    ; VEGF  : sigmoid Imax drug effect relationship
 EFF2    =   IMAX2*AUC**HILL2/(IC502**HILL2+AUC**HILL2) ; VEGFR2: sigmoid Imax drug effect relationship
 EFF3    =   IMAX3*AUC/(IC503+AUC)                      ; VEGFR3: Imax drug effect relationship
 EFFS    =   IMAXS*AUC/(IC50S+AUC)                      ; SKIT  : Imax drug effect relationship


 DADT(1) =   KIN-KOUT*(1-EFF)*A(1)                      ; VEGF  : Inhibition of Kout
 DADT(2) =   KIN2*(1-EFF2)-KOUT2*A(2)                   ; VEGFR2: Inhibition of Kin
 DADT(3) =   KIN3*(1-EFF3)-KOUT3*A(3)                   ; VEGFR3: Inhibition of Kin
 DADT(4) =   KINS*(1-EFFS)-KOUTS*A(4)                   ; SKIT  : Inhibition of Kin


$ERROR

;----- stratification for VPCs ----------
STRT=0
IF(FLAG.EQ.5.AND.PLA.EQ.0)  STRT=1                      ; VEGF, treatment
IF(FLAG.EQ.5.AND.PLA.EQ.1)  STRT=2                      ; VEGF, placebo
IF(FLAG.EQ.6)               STRT=3                      ; sVEGFR-2
IF(FLAG.EQ.7)               STRT=4                      ; sVEGFR-3
IF(FLAG.EQ.8.AND.PLA.EQ.0)  STRT=5                      ; sKIT, treatment
IF(FLAG.EQ.8.AND.PLA.EQ.1)  STRT=6                      ; sKIT, placebo


  DEL=0                                                 ; to avoid numerical problems
  IF(A(1).EQ.0)DEL=1
  IF(A(2).EQ.0)DEL=1
  IF(A(3).EQ.0)DEL=1
  IF(A(4).EQ.0)DEL=1

  IF(FLAG.EQ.5) THEN                                    ; VEGF
  IPRED  =  LOG(A(1)+DEL)                               ; modeling performed on log transformed VEGF
  W      =  THETA(14)
  Y      =  IPRED+W*EPS(1)
  ENDIF

  IF(FLAG.EQ.6) THEN                                    ; VEGFR2
  IPRED  =  LOG(A(2)+DEL)                               ; modeling performed on log transformed VEGFR2
  W      =  SQRT(THETA(15)**2+(THETA(16)/A(2))**2)
  Y      =  IPRED+W*EPS(1)
  END IF


  IF(FLAG.EQ.7) THEN                                    ; VEGFR3
  IPRED  =  LOG(A(3)+DEL)                               ; modeling performed on log transformed VEGFR3
  W      =  THETA(17)
  Y      =  IPRED+W*EPS(1)
  ENDIF

  IF(FLAG.EQ.8) THEN                                    ; SKIT
  IPRED  =  LOG(A(4)+DEL)                               ; modeling performed on log transformed SKIT
  W      =  THETA(18)
  Y      =  IPRED+W*EPS(1)
  ENDIF


  IF(W.EQ.0)W=1
  IRES   =  DV -  IPRED
  IWRES  =  IRES/W

;-----INITIAL ESTIMATES------------------------

;-----VEGF-------------------------------------
$THETA   (0,59.7)                          ; 1  BM0
$THETA   (0,91)                            ; 2  MRT
$THETA   1 FIX                             ; 3  IMAX
$THETA   (0,1)                             ; 4  IC50
$THETA   (0,3.31)                          ; 5  HILL
$THETA   (-0.06,0.035)                     ; 6  DP SLOPE
;-----VEGFR2-----------------------------------
$THETA   (0,8670)                          ; 7  BM02
$THETA   (0,554)                           ; 8  MRT2
$THETA   (0,1.54)                          ; 9  HILL2
;-----VEGR3------------------------------------
$THETA   (0,63900)                         ; 10 BM03
$THETA   (0,401)                           ; 11 MRT3
;-----SKIT-------------------------------------
$THETA   (0,39200)                         ; 12 BM0S
$THETA   (0,2430)                          ; 13 MRTS

;-----RES--------------------------------------

$THETA   (0,0.445)                         ; 14 RES ERROR PROP
$THETA   (0,0.12)                          ; 15 RES VEGFR2
$THETA   (0,583)                           ; 16 RES VEGFR2
$THETA   (0,0.22)                          ; 17 RES VEGFR3
$THETA   (0,0.224)                         ; 18 RES sKIT
;-----BM0-------------------------------------
$OMEGA   0.252                             ; 1  ETA BM0
$OMEGA   0.0369                            ; 2  ETA BM02
$OMEGA   0.186                             ; 3  ETA BM03
$OMEGA   0.254                             ; 4  ETA BM0S
;-----MRT--------------------------------------
$OMEGA   0.0600                            ; 5  ETA MRT sVEGFR-2,3
$OMEGA   0.0753                            ; 6  ETA MRT sKIT

;-----IC50-------------------------------------
$OMEGA BLOCK(4)   0.253                    ; 7  ETA IC50
                  0.198 0.189              ; 8  ETA IC502
                  0.238 0.252 0.398        ; 9  ETA IC503
                  0.218 0.297 0.936 5.77   ; 10 ETA IC50S

$OMEGA   2.95                              ; 11 ETA DP SLOPE VEGF
$OMEGA   3.01                              ; 12 ETA DP SLOPE SKIT

$SIGMA   1 FIX

$EST SIGDIGITS=4 MAXEVALS=0 POSTHOC PRINT=1 MSFO=MSF1 METH=1 INTER        ; FOCE with interaction
;$COV
$TABLE ID CYCL EVID TIME PLA FLAG IPRED IWRES NOPRINT ONEHEADER FILE=sdtab1
$TABLE ID FLAG BM0 MRT IC50 DPSLO BM02 MRT2 IC502 BM03 MRT3 IC503 BM0S MRTS IC50S DPSLOS NOPRINT ONEHEADER FILE=patab1
$TABLE ID TIME PLA FLAG  NOPRINT ONEHEADER FILE=catab1