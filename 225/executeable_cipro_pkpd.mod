
$PROBLEM Khan et al, JAC 2015

$INPUT STR ID MIC xMIC CAB TIME
       TUBE CMT EVID DV AMT FLG2 FLG1 BASE L2 BOBS

$DATA ciprotest.csv IGNORE=@  IGNORE=(BOBS.EQ.1)   ;BOBS=baseline observaions. Here dropped, using B2 for baseline.

$SUBS  ADVAN13 TOL5

$MODEL
NCOMP=6
COMPARTMENT=(S)       ;Compartment S1
COMPARTMENT=(R)       ;Compartment R1
COMPARTMENT=(SPE)     ;Compartment S2
COMPARTMENT=(NP)      ;Compartment NC1
COMPARTMENT=(RPE)     ;Compartment R2
COMPARTMENT=(NPPE)    ;Compartment NC2

$PK
" FIRST
"  COMMON /PRCOMG/ IDUM1,IDUM2,IMAX,IDUM4,IDUM5
"  INTEGER IDUM1,IDUM2,IMAX,IDUM4,IDUM5
"  IMAX=10000000


;;-- growth rate

KGS   = THETA(1)

;;-- death rate
KK    = 0.179

;;-- drug EMAX

EMAX  = THETA(2)

;;-- drug EC50
IF(STR.EQ.347)   EC50  = THETA(3)
IF(STR.EQ.202)   EC50  = THETA(4)
IF(STR.EQ.378)   EC50  = THETA(5)
IF(STR.EQ.534)   EC50  = THETA(6)
IF(STR.EQ.625)   EC50  = THETA(7)
IF(STR.EQ.693)   EC50  = THETA(8)
IF(STR.EQ.707)   EC50  = THETA(9)


;;-- hill factor drug effect

GAM   = THETA(10)

;; -- proportionality constant

PC  = THETA(11)*0.0000001     ; proportionality constant

;; -- baseline, B2 method

SBASE=EXP(BASE)*EXP(ETA(1)*(SIGMA(1)**0.5))

;;--growth rate 2

KGS2= THETA(12)

;;-- drug EC502
EC502= THETA(13)

;;--pre-existing resistant at start
MUT= THETA(14)

;;--compartment initialization

A_0(1)= SBASE*(1-MUT*0.000001)           
A_0(2)= 0                                
A_0(3)= MUT*0.000001*SBASE
A_0(4)= 0
A_0(5)= 0
A_0(6)= 0


;;-- k non-colony forming compartements

KSNC1    =  THETA(15) * (CAB/EC50)**THETA(17)/((CAB/EC50)**THETA(17)+THETA(18)**THETA(17))
KSNC2   =  THETA(15) * (CAB/EC502)**THETA(17)/((CAB/EC502)**THETA(17)+THETA(18)**THETA(17))

KNCS1    =  THETA(16)* EC50/(CAB+0.0000000001)   ; avoids division by 0
KNCS2   =  THETA(16)* EC502/(CAB+0.0000000001)


;;-- MTIME - shuts off KNC1
MTIME(1) = 0
MTIME(2) = MTIME(1)+THETA(19)


$DES

;;-- see MTIME

FLAG=MPAST(1)-MPAST(2)

;;-- Conversion rate between active/resting cell dependent on cell number

SR=PC*(A(1)+A(2)+A(3)+A(4)+A(5)+A(6))
SR2=SR
RS=0
RS2=0


;EMAX equation for drug effect

DRUGS=0
DRUGS2=0
IF(CAB.GT.0.00000000001)THEN
 DRUGS = EMAX*(CAB)**GAM/(EC50**GAM+(CAB)**GAM)
 DRUGS2 = EMAX*(CAB)**GAM/(EC502**GAM+(CAB)**GAM)
ENDIF

;;-- Dif.eq.
DADT(1)= KGS*(A(1))-(KK+DRUGS)*(A(1)) - SR*(A(1)) + RS*(A(2)) +KNCS1*A(4) -KSNC1*A(1)*FLAG          ; S
DADT(2)=-KK*(A(2)) + SR*(A(1)) - RS*(A(2))                                                          ; R
DADT(3)= KGS2*(A(3))-(KK+DRUGS2)*(A(3)) - SR2*(A(3)) + RS2*(A(5))+KNCS2*A(6) -KSNC2*A(3)*FLAG     ; S2
DADT(4)= KSNC1*A(1)*FLAG - KNCS1*A(4)-(KK+DRUGS)*(A(4))                                             ; NC1
DADT(5)=-KK*(A(5)) + SR2*(A(3))                                                                     ; R2
DADT(6)= KSNC2*A(3)*FLAG - KNCS2*A(6)-(KK+DRUGS2)*(A(6))                                          ; NC2

$ERROR
SIG=SQRT(SIGMA(1)+SIGMA(2))
LLOQ=LOG(10)  ;LOQ=10 (LOG=LN in NONMEM)

A1=A(1)
A2=A(2)
A3=A(3)
A4=A(4)
A5=A(5)
A6=A(6)

ATOT=A1+A2+A3+A5       ;NC not counted
DEL=0
IPRED=LOG(ATOT+1E-8)

DUM = (LLOQ-IPRED)/SIG
DUM2= PHI(DUM)

IRES  = DV-IPRED
W=SQRT(SIGMA(1)+SIGMA(2))
IWRES = IRES/W

SAM1=0
SAM2=0
SAM3=0
SAM4=0
IF(FLG2.EQ.1) SAM1=1
IF(FLG2.EQ.2) SAM2=1
IF(FLG2.EQ.3) SAM3=1
IF(FLG2.EQ.4) SAM4=1

IF(FLG1.EQ.0) THEN     ;DV>LOQ
F_FLAG=0
Y=IPRED+EPS(1)+SAM1*EPS(2)+SAM2*EPS(3)+SAM3*EPS(4)+SAM4*EPS(5) 
ENDIF

IF(FLG1.EQ.1) THEN     ;DV<LOQ
F_FLAG=1
Y=DUM2
ENDIF


$THETA  (1,1.70)        ; 1. kgrowth1
$THETA  (0,5.24)        ; 2. EMAX
$THETA  (0,0.037) FIX   ; 3. Ec50 347
$THETA  (0,0.057)       ; 4. Ec50 202  only LM202 provided in ciprotest.csv dataset                                                           '
$THETA  (0,0.65)  FIX   ; 5. Ec50 378
$THETA  (0,0.30)  FIX   ; 6. Ec50 534
$THETA  (0,1.0)   FIX   ; 7. Ec50 625
$THETA  (0,31)    FIX   ; 8. Ec50 693
$THETA  (0,92)    FIX   ; 9. Ec50 707
$THETA  (0.5,1.98)   ; 10. GAM
$THETA  (0,0.0186)   ; 11. PC
$THETA  (0.18,0.344) ; 12. kgrowth2
$THETA  (0,1.25)     ; 13. EC502
$THETA  (0,0.81)     ; 14. start conc pre-existing resistant
$THETA  (0,5.83)     ; 15. kSNc,max
$THETA  (0,0.17)     ; 16. sfNcS
$THETA  (0,20) FIX   ; 17. Hill factor Nc
$THETA  (0,0.24)     ; 18. tr50
$THETA  (2,5.3158)   ; 19. MTIME

$OMEGA  1  FIX  ; B2 method OMEGA 1 FIX

$SIGMA  2.41597
$SIGMA  BLOCK(1) 0.128122
$SIGMA  BLOCK(1) SAME
$SIGMA  BLOCK(1) SAME
$SIGMA  BLOCK(1) SAME

$ESTIMATION METHOD=1 LAPLACIAN NUMERICAL NSIG=5   MAXEVAL=9999 POSTHOC  PRINT=1 NOABORT  INTERACTION
;$COVARIANCE  UNCONDITIONAL
$TABLE ID TIME CAB IPRED IWRES PRED NOPRINT ONEHEADER FILE=sdtab1
$TABLE SBASE DV BASE DRUGS ATOT A1 A2 A3 A4 CAB ID TIME KGS SR RS KK EMAX                        
 KGS2 EC502 PRED GAM FLAG IPRED STR DRUGS2 NOPRINT ONEHEADER FILE=patab1
