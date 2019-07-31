$PROBLEM The Multistate Tuberculosis Pharmacometric Model
;The Multistate Tuberculosis Pharmacometric Model
;Oskar Clewe
;Uppsala University
;October 2015
;M.tuberculosis H37Rv (St George´s Uni)
;NM7.3
;Log-transformed both sides
;DV is ln(CFU/mL)
;Time in days
$INPUT      TIME ID NDV DV EVID MDV AMT
$DATA      Data.csv IGNORE=@
$SUBROUTINE ADVAN13 TOL=9
$MODEL      NCOMP=3 COMP=(FBUGS) COMP=(SBUGS) COMP=(NBUGS)
$PK
TVKG=THETA(1)            ; Growth rate of F bacteria

KFSLIN=THETA(2)/100    ; Rate parameter, F -> S, Linear time dependent
KFN=THETA(3)/1000000   ; Rate parameter, F -> N
KSF=THETA(4)/10        ; Rate parameter, S -> F
KSN=THETA(5)           ; Rate parameter, S -> N
KNS=THETA(6)/100       ; Rate parameter, N -> S

TVBMAX=THETA(7)*1000000  ; System carrying capacity (CFU/ml)

TVF0=THETA(8)            ; Initial F bacterial number (CFU/ml)
TVS0=THETA(9)            ; Initial S bacterial number (CFU/ml)

KG=TVKG
BMAX=TVBMAX
F0=TVF0*EXP(ETA(1))      ; IIV on initial F bacterial number
S0=TVS0

A_0(1)=F0                ; Initial F bacterial number with IIV
A_0(2)=S0                ; Initial S bacterial number
A_0(3)=0.00001           ; Initial N bacterial number

$DES       
GROWTHFUNC=KG*LOG(BMAX/(A(1)+A(2)+A(3))) ; Gompertz growth function
IF(GROWTHFUNC.LT.0) GROWTHFUNC=0         ; Keep GROWTHFUNC from turning negative

KFS=KFSLIN*T                             ; Linear time dependendent transfer, F -> S

DADT(1)=A(1)*GROWTHFUNC+KSF*A(2)-KFS*A(1)-KFN*A(1) ;F
DADT(2)=KFS*A(1)+KNS*A(3)-KSN*A(2)-KSF*A(2)        ;S
DADT(3)=KSN*A(2)+KFN*A(1)-KNS*A(3)                 ;N

$ERROR       
FBUGS=A(1)             ; F
SBUGS=A(2)             ; S
NBUGS=A(3)             ; N
TOTBUGS=A(1)+A(2)+A(3) ; F+S+N

IPRED=LOG(A(1)+A(2))
IRES=DV-IPRED
ADD=SQRT(SIGMA(1))
SD=SQRT((ADD)**2)      ; SD for additive residual error on log scale
IWRES=IRES/SD
Y=IPRED+EPS(1)

$THETA  (0,0.206361)   ; 1 kG
$THETA  (0,0.1657)     ; 2 kFSLIN (/100)
$THETA  (0,0.9)        ; 3 kFN (/1000000)
$THETA  (0,0.14478)    ; 4 kSF (/10)
$THETA  (0,0.185568)   ; 5 kSN
$THETA  (0,0.1227)     ; 6 kNS (/100)
$THETA  (0,241.6170)   ; 7 Bmax (*1000000)
$THETA  (0,4.109880)   ; 8 F0
$THETA  (0,9770.730)   ; 9 S0
$OMEGA  22.37250       ; variance for ETA(1) on F0
$SIGMA  0.400262       ; variance for add residual error on logscale
$ESTIMATION METHOD=1 MAXEVAL=9999 NSIG=3 SIGL=9
$COVARIANCE PRINT=E
$TABLE      ID TIME IPRED ADD IRES IWRES CWRES DV NDV FBUGS SBUGS
            NBUGS TOTBUGS EVID ONEHEADER NOPRINT FILE=sdtabMTP
$TABLE      ID TIME GROWTHFUNC KG KFN KFS KFSLIN KSF KSN KNS
            BMAX F0 S0 ETA(1) ONEHEADER NOPRINT FILE=patabMTP
$TABLE      ID TIME ONEHEADER NOPRINT FILE=cotabMTP
$TABLE      ID TIME ONEHEADER NOPRINT FILE=catabMTP

