$PROBLEM    RIF final model with simulated data
$INPUT      ID OCC TIME DV MDV EVID DOSE AMT FDC
$DATA       Simulated_TB_Rifampicin_PK_Wilkins_2008.csv IGNORE=@
$SUBROUTINE ADVAN6 TRANS1 TOL=4
$MODEL      NCOMP=2 COMP=(DEPOT,DEFDOSE) COMP=(CENTRAL)
$PK  

"FIRST
"        COMMON/PRCOMG/  IDUM1,IDUM2,IMAX,IDUM4,IDUM5
"INTEGER IDUM1,IDUM2,IMAX,IDUM4,IDUM5
"IMAX=10000000

;; covariates
IF(FDC.EQ.1) MTTFDC = 0  
IF(FDC.EQ.0) MTTFDC = THETA(8)

MTTCOV = (1 + MTTFDC)

IF(FDC.EQ.1) CLLOC = 0 
IF(FDC.EQ.0) CLLOC = THETA(9)

CLCOV = (1 + CLLOC)

;; set up dosing and time after dose
IF(AMT.GT.0) PD   = AMT
IF(AMT.GT.0) TDOS = TIME
TAD = TIME-TDOS

IOVCL = ETA(5)
IOVMT = ETA(11)
IF (OCC.EQ.2) THEN
  IOVCL = ETA(6)
  IOVMT = ETA(12)
ENDIF
IF (OCC.EQ.3) THEN
  IOVCL = ETA(7)
  IOVMT = ETA(13)
ENDIF
IF (OCC.EQ.4) THEN 
  IOVCL = ETA(8)
  IOVMT = ETA(14)
ENDIF
IF (OCC.EQ.5) THEN 
  IOVCL = ETA(9)
  IOVMT = ETA(15)
ENDIF
IF (OCC.EQ.6) THEN
  IOVCL = ETA(10)
  IOVMT = ETA(16)
ENDIF

IIVCL  = ETA(1)
IIVMT  = ETA(4)

TVCL   = THETA(1)*CLCOV              ; CL
CL     = TVCL*EXP(IOVCL + IIVCL)     ; CL

TVV    = THETA(2)                    ; V
V      = TVV*EXP(ETA(2))             ; V

K     = CL/V                         ; K

TVKA  = THETA(3)                     ; KA
KA    = TVKA * EXP(ETA(3))           ; KA

S2    = V
F1  = 0                            ; important!

TVMTT = THETA(6)*MTTCOV              ; MTT
MTT   = TVMTT*EXP(IOVMT + IIVMT)     ; MTT

TVNN  = THETA(7)                     ; NN
NN    = TVNN*EXP(ETA(17))            ; NN

KTR   = (NN+1)/MTT
L     = LOG(2.5066)+(NN+.5)*LOG(NN)-NN+LOG(1+1/(12*NN))

$DES  

X = 0.00001

IF(T.GT.TDOS)THEN
  DADT(1)=EXP(LOG(PD+X)+LOG(KTR+X)+NN*LOG(KTR*(T-TDOS)+X)-KTR*(T-TDOS)-L)-KA*A(1)
ELSE
  DADT(1)=EXP(LOG(PD+X)+LOG(KTR+X)+NN*LOG(KTR*T+X)-KTR*T-L)-KA*A(1)
ENDIF
DADT(2)=KA*A(1)-K*A(2)

$ERROR   (ONLY OBSERVATIONS)

AA1   = A(1)
AA2   = A(2)

CP    = A(2)/V

IPRED = F
W     = SQRT(THETA(4)**2+THETA(5)**2*F*F)
IRES  = DV-IPRED
IWRES = IRES/W
Y     = IPRED+W*EPS(1)

$THETA  (0,20.4937) ; 1 CL
$THETA  (0,55.8433,100) ; 2 V
$THETA  (0,0.878858) ; 3 KA
$THETA  (0,0.0803316) ; 4 ADD err
$THETA  (0,0.255282) ; 5 CCV err
$THETA  (0,0.382318) ; 6 MTT
$THETA  (1,4.42696,80) ; 7 NN
$THETA  (-1,0.935562,5) ; 8 FDC on MTT
$THETA  (-1,0.155522,5) ; 9 LOC on CL
$OMEGA  BLOCK(2)
 0.237494  ; 1 IIV CL/F
 0.171532 0.135696  ;  2 IIV V/F
$OMEGA  0.454922  ;   3 IIV KA
$OMEGA  0.333651  ;  4 IIV MTT
$OMEGA  BLOCK(1)
 0.048617  ;   5 IOV CL
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1)
 0.253377  ; 11 IOV MTT
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  0.959875  ;      17 NN
$SIGMA  1.000000  FIX
$ESTIMATION METHOD=1 INTER PRINT=3 MAXEVAL=9999 SIGDIG=3
            MSFO=TB_Rifampicin_PK_Wilkins_2008_simulated.msf
$COVARIANCE PRINT=E MATRIX=S
$TABLE      ID TIME IPRED IWRES AA1 AA2 CP NOPRINT ONEHEADER
            FILE=sdtab.simulated_TB_Rifampicin_PK_Wilkins_2008
$TABLE      ID CL V KA K MTT NN ETA1 ETA2 ETA3 ETA4 ETA5 ETA11 ETA17
            NOPRINT ONEHEADER FILE=patab.simulated_TB_Rifampicin_PK_Wilkins_2008
$TABLE      ID DOSE NOPRINT ONEHEADER
            FILE=cotab.simulated_TB_Rifampicin_PK_Wilkins_2008
$TABLE      ID FDC NOPRINT ONEHEADER
            FILE=catab.simulated_TB_Rifampicin_PK_Wilkins_2008

