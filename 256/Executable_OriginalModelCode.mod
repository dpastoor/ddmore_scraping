$PROBLEM Phenobarbital PK in newborns 
$INPUT   ID 
         DOSE=AMT ;in mg
         CONC=DV  ;in mg/L
         WEIGHT   ;in kg
         CMT      ;1 oral depot, 2 central
         TIME     ;in hours
         RATE     ;in mg/h
         BWEIGHT  ;birthweight in kg   
         AGE      ;postnatal age in days
         MDV
$DATA SimulatedPhenobarbitalNewbornsPK.csv IGNORE=@
$SUBROUTINE ADVAN2 TRANS2 
$PK 

;;; VWEIGHT-DEFINITION START 
VWEIGHT = ( 1 + THETA(6)*(WEIGHT - 2.70)) 
;;; VWEIGHT-DEFINITION END 

;;; V-RELATION START 
VCOV=VWEIGHT 
;;; V-RELATION END 


;;; CLBW-DEFINITION START 
CLBW = ( 1 + THETA(5)*(BWEIGHT - 2.59)) 
;;; CLBW-DEFINITION END 

;;; CLAGE-DEFINITION START 
CLAGE = ( 1 + THETA(4)*(AGE - 4.50)) 
;;; CLAGE-DEFINITION END 

;;; CL-RELATION START 
CLCOV=CLAGE*CLBW 
;;; CL-RELATION END 


TVCL = THETA(1) * CLCOV
TVV  = THETA(2) * VCOV 
CL   = TVCL*EXP(ETA(1)) 
V    = TVV*EXP(ETA(2))
F1   = THETA(8)
KA   = THETA(7)
S2   = V
K    = CL/V 

$ERROR  
IPRED = F 
W = SQRT(THETA(3) * IPRED**2) 
IWRES = (DV-IPRED) / W 
Y= IPRED + EPS(1) * W

$THETA
(0, 0.00909) ;CL 
(0, 2.38) ;V 
(0, 0.0258) ; prop error 
(-0.01, 0.0533,0.2) ; CLAGE 
(-0.205, 0.369,0.5) ; CLBW 
(-0.555, 0.309,0.444) ; VWEIGHT 
(50) FIX ; KA
(0, 0.594,1) ; F

$OMEGA 0.0898
$OMEGA 0.0504

$SIGMA 1 FIX 

$ESTIMATION METHOD=1 INTERACTION MAXEVALS=9990 POSTHOC
$COVARIANCE 
$TABLE ID TIME PRED IPRED RES WRES IWRES CWRES K KA CL V ETA1 ETA2
ONEHEADER NOPRINT NOAPPEND FILE=sdtab524

