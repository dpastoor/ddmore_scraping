;; 1. Based on: run126g
;; 2. Description: For DDMORE model repository
;; x1. Author: user
$PROBLEM   
$INPUT      ID TIME STUDY EVID CMT AMT RATE DV T2DM BASI INSU BW
$DATA      Simulated_ddmoremockdata2.txt IGNORE=@
$SUBROUTINE ADVAN13 TRANS1 TOL=6
$MODEL      COMP=STOMACH ;1 Stomach compartment
COMP=INTESTI ;2 Intestine compartment
COMP=CENTP ;3 Central compartment
COMP=PERIPHE ;4 Peripheral compartment
COMP=SG ;5 Glucose stomach
COMP=DG ;6 GLucose duodenum
COMP=CENTG ;7 Central glucose
COMP=PERIG ;8 Peripheral glucose
COMP=EFF_G ;9 Effect compartment
COMP=EFF_I ;10 Effect compartment
COMP=JEJUN ;11 Jejunum compartment
COMP=ILEUM ;12 Ileum compartment
COMP=KUML ;13 Cumulative first pass loss paracetamol
COMP=GLP1 ;14 GLP1
COMP=GIP ;15 GIP
$PK 
;---- Paracetamol parameters ---
CL     = THETA(1)*EXP(ETA(1)) ; Paracetamol clearance
V1     = THETA(2)*EXP(ETA(2)) ; Central volume of distribution
Q2     = THETA(3)             ; Intercompartmental clearance
V2     = THETA(4)             ; Peripheral volume of distribution
KA     = THETA(5)             ; Absorption rate constant

IF(STUDY.EQ.1) APAPBL = THETA(6)*EXP(ETA(3))       ; Baseline noise
IF(STUDY.EQ.2) APAPBL = THETA(6)*EXP(ETA(4))       ; Baseline noise
IF(STUDY.EQ.3) APAPBL = THETA(6)*EXP(ETA(5))       ; Baseline noise

K12A   = Q2/V1
K21A   = Q2/V2
KE     = CL/V1

;---- Glucose parameters -------
GSSH   = THETA(7)*EXP(ETA(6))   ; Glucose baseline healthy
GSSD   = THETA(8)*EXP(ETA(7))   ; Glucose baseline T2DM

VG     = THETA(9)*BW/70         ;Central volume for glucose
CLGH   = THETA(10)*EXP(ETA(8)) ;Glucose clearance
CLGIH  = THETA(11)*EXP(ETA(9)) ;Insulin dependant glucose clearance
Q      = THETA(12)              ;Intercompartment glucose clearance
VP     = THETA(13)              ;Periferal volume for glucose
KGE1   = THETA(14)
KIE    = THETA(15)
CLGD   = THETA(16)*EXP(ETA(10))
CLGID  = THETA(17)*EXP(ETA(11))

FPGH  = THETA(18)
FPGD  = THETA(19)
RAMAXD = THETA(20)  ; Maximum rate of absorption glucose
RAMAXJ = THETA(21)
RAMAXI = THETA(22)
KMG     = THETA(23)  ; 50 % of absorption rate

;-------------------------------
IF(T2DM.EQ.0) GSS = GSSH ; Healthy have healthy glucose baseline
IF(T2DM.EQ.1) GSS = GSSD ; T2DM have T2DM glucose baseline
IF(T2DM.EQ.0) CLG = CLGH ; Healthy
IF(T2DM.EQ.1) CLG = CLGD ; T2DM
IF(T2DM.EQ.0) CLGI = CLGIH ; Healthy
IF(T2DM.EQ.1) CLGI = CLGID ; T2DM
IF(T2DM.EQ.0) GPRG = -2.79 ; Healthy
IF(T2DM.EQ.1) GPRG = 0     ; T2DM
IF(T2DM.EQ.0) FPG = FPGH ; Healthy
IF(T2DM.EQ.1) FPG = FPGD ; T2DM

ISS  = BASI/6.945 ; Obeserved baseline in mU/L

K12G = Q/VG
K21G = Q/VP
KG   = CLG/VG
KGI  = CLGI/VG

GPRO = GSS*(KG+KGI*ISS)*VG*180/1000     ;Baseline glucose production in amount mMol => g, since dose is in g

;---- Gastric emptying ---------------
KW     = THETA(24)*EXP(ETA(12))  ; "Water" emptying
GLUMAX   = THETA(25)              ; Maximum inhibition of glucose dose
IGD50  = THETA(26)*EXP(ETA(13)) ; 50% of glucose inhibition
GAM      = THETA(27)
IF(STUDY.EQ.1) T50 = THETA(28)*EXP(ETA(14))
IF(STUDY.EQ.2) T50 = THETA(28)*EXP(ETA(15))
IF(STUDY.EQ.3) T50 = THETA(28)*EXP(ETA(16))

;----- Saturable first pass effect ---
VMAX   = THETA(29)
KMAC = THETA(30)

;----- Incretin parameters -----------
EMGLP1D  = THETA(31)
EMGLP1I  = THETA(32)
E50GLP1I = THETA(33)

EMGIP   = THETA(34)
E50GIP  = THETA(35)

BGLP1 = THETA(36)*EXP(ETA(17))
BGIP  = THETA(37)*EXP(ETA(18))

KEGIP = LOG(2)/7
KEGLP1 = LOG(2)/1.5

;----- GI parameters -----------------
TT   = 240 ; min transit time through small intestine

FLD  = 0.08 ; Fractional length duodenum
FLJ  = 0.37 ; Fractional length jejunum
FLI  = 0.55 ; Fractional length ileum

DTT  = TT*FLD ; Duodenum transit time
KDJ  = 1/DTT ; K duodenum jejunum

JTT  = TT*FLJ ; Jejunum transit time
KJI  = 1/JTT ; K jejunum ileum

ITT  = TT*FLI ; Ileum transit time
KIC  = 1/ITT ; K ileum colon

;-----------------------------------------------------
A_0(1) = 0                         ; Stomach paracetamol
A_0(2) = 0                         ; Duodenum paracetamol
A_0(3) = 0                         ; Central paracetamol
A_0(4) = 0                         ; Peripheral paracetamol
A_0(5) = 0                         ; Stomach glucose
A_0(6) = 0                         ; Duodenum glucose
A_0(7) = GSS*VG*180/1000           ; Central glucose
A_0(8) = K12G/K21G*GSS*VG*180/1000 ; Peripheral glucose
A_0(9) = GSS                       ; Glucose on glucose production
A_0(10)= ISS                       ; Insulin on glucose elimination
A_0(11) = 0                        ; Jejunum glucose
A_0(12) = 0                        ; Ileum glucose
A_0(13) = 0                        ; CFPL
A_0(14) = BGLP1
A_0(15) = BGIP

$DES
GLURAT = A(9)/GSS
GGPR   = (GLURAT+0.0001)**GPRG   ; Glucose on glucose production
GPRT   = GPRO*GGPR               ; Glucose production   a.k.a GPRD

RAD    = RAMAXD*A(6)/(KMG+A(6))
RAJ    = RAMAXJ*A(11)/(KMG+A(11))
RAI    = RAMAXI*A(12)/(KMG+A(12))
RA     = RAD + RAJ + RAI

DUOGLU = A(6)
IF(DUOGLU.LT.0) DUOGLU = 0
GLUINHI = GLUMAX*(DUOGLU+0.00001)**GAM/(IGD50**GAM+(DUOGLU+0.00001)**GAM)
KS = KW*(1-GLUINHI)
IF(KS.LT.0) KS = 0

LAG = 1/(1+EXP(-10*(T-T50)))
FPL = VMAX*A(2)/(KMAC+A(2))

DADT(1)  = -KS*A(1)*LAG                             ;1 Stomach
DADT(2)  = KS*A(1)*LAG -KA*A(2)-FPL                     ;2 Intestine
DADT(3)  = KA*A(2)- (KE+K12A)*A(3) + K21A*A(4);3 Central
DADT(4)  = K12A*A(3) - K21A*A(4)                    ;4 Peripheral

DADT(5)  = -KS*A(5)*LAG                             ;5 Stomach compartment (g)
DADT(6)  = KS*A(5)*LAG -RAD -KDJ*A(6)               ;6 Duodenum compartment (g)
DADT(7)  = RA*FPG + GPRT +K21G*A(8) -(K12G+KG)*A(7)- KGI*A(7)*A(10) ;7 Central glucose, amount (g)
DADT(8)  = K12G*A(7) -K21G*A(8)                     ;8 Periferal glucose (g)
DADT(9)  = KGE1*(A(7)/VG/180*1000 - A(9))           ;9 Effect compartment, Glucose on glucose production g => mM
DADT(10) = KIE*(INSU/6.945 - A(10))                 ;10 Effect compartment, Insulin on glucose elimination, pM => µU/ml
DADT(11) = KDJ*A(6)  -KJI*A(11) - RAJ               ;11 Jejunum
DADT(12) = KJI*A(11) - RAI               ;12 Ileum
DADT(13) = FPL                          ;13 Cumulative amount lost

;----- GIP & GLP1 ------------------------------
DADT(14) = -KEGLP1*(A(14)-BGLP1)+ BGLP1*KEGLP1*(EMGLP1D*A(6)+EMGLP1I*A(12)/(E50GLP1I+A(12)))   ; 14 GLP1
DADT(15) = -KEGIP *(A(15)-BGIP )+ BGIP *KEGIP *              EMGIP*A(6)/(A(6)+E50GIP)          ; 15 GIP
;-----------------------------------------------------------

$ERROR 
STOM  = A(1)
INTE  = A(2)
CAC   = A(3)/151.17/V1*1000 ; Convert amount (mg) to concentration (µM)
PAC   = A(4)
STOMG = A(5)
DUODG = A(6)
CGLU  = A(7)*1000/(180*VG) ; convert g to mM
CPER  = A(8)*1000/(180*VP) ; convert g to mM
GGP   = A(9)
IGE   = A(10)
JEJUG = A(11)
ILEUG = A(12)
CFPL  = A(13)
CGLP1 = A(14)
CGIP  = A(15)

IF(CMT.EQ.3) IPRED = LOG(CAC+APAPBL+0.00001)
IF(CMT.EQ.3) W = SQRT(SIGMA(1))
IF(CMT.EQ.3) Y = IPRED + EPS(1)


IF(CMT.EQ.7) IPRED = LOG(CGLU+0.00001)
IF(CMT.EQ.7) W = SQRT(SIGMA(2))
IF(CMT.EQ.7) Y = IPRED + EPS(2)

IF(CMT.EQ.14) IPRED = LOG(CGLP1+0.00001)
IF(CMT.EQ.14) W = SQRT(SIGMA(3))
IF(CMT.EQ.14) Y = IPRED + EPS(3)

IF(CMT.EQ.15) IPRED = LOG(CGIP+0.00001)
IF(CMT.EQ.15) W = SQRT(SIGMA(4))
IF(CMT.EQ.15) Y = IPRED + EPS(4)

IRES  = DV-IPRED
IF(W.EQ.0) W = 1
IWRES = IRES/W

;--- Initial estimates ---
$THETA
(0, 0.344) ; 1 CL (L/min)
(0, 27) ; 2 V1 (L)
(0, 0.675) ; 3 CLQ(L/min)
(0, 22) ; 4 V2 (L)
(0.14) FIX ; 5 KA (/min) 5 min
(0, 6.54) ; 6 APAPBL (µM)
(0, 5.27) ; 7 GSSH   (mM)
(0, 7.48) ; 8 GSSD   (mM)
(9.33) FIX ; 9 VG (L)
(0.0894) FIX ; 10 CLGH (L/min)
(0, 0.00663,0.0083) ; 11 CLGIH (L/min/(µU/mL)
(0.442) FIX ; 12 Q (L/min)
(8.56) FIX ; 13 VP (L)
(0.0573) FIX ; 14 KGE1
(0.0213) FIX ; 15 KIE
(0.0287) FIX ; 16 CLGD (L/min)
(0, 0.00547,0.0083) ; 17 CLGID (L/min/(µU/mL)
(0, 0.909,1) ; 18 FPGH
(1) FIX ; 19 FPGD
(0, 0.576) ; 20 RAMAXD rate
(0, 2.06) ; 21 RAMAXJ rate
(0, 1.33) ; 22 RAMAXI rate
(0, 6.32) ; 23 KMG (g)
(0.14) FIX ; 24 KW   (/min)
(1) FIX ; 25 GLUMAX
(0, 7.42) ; 26 IGD50
(0, 14) ; 27 GAM
(15, 20,30) ; 28 T50
(0, 17) ; 29 VMAX
(0, 168) ; 30 KMAC
(0, 0.0279,300) ; 31 EMGLP1D
(0, 1.98,300) ; 32 EMGLP1I
(0, 0.0525) ; 33 E50GLPI
(0, 14.3,300) ; 34 EMGIP
(0, 2.57,300) ; 35 E50GIP
(0, 17.1) ; 36 BGLP1
(0, 8.03) ; 37 BGIP

$OMEGA 0.109 ;   1 IIV CL
$OMEGA 0.143 ;   2 IIV V1
$OMEGA  BLOCK(1) 0.0338 ;   3 APAPBL
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
;Glucose
$OMEGA 0.00177 ;     6 GSSH
$OMEGA 0.00762 ;     7 GSSD
$OMEGA 0.348 FIX  ;     8 CLGH
$OMEGA 0.226 ;    9 CLGIH
$OMEGA 0.348 FIX  ;    10 CLGD
$OMEGA 0.209 ;   11 CLGID
;Gastric emptying
$OMEGA 0.25 FIX  ;    12 KW
$OMEGA 0.114 ; 13 IGD50
$OMEGA  BLOCK(1) 0.0538 ;     14 LAG
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
;Incretins
$OMEGA 0.0649 ;   17 BGLP1
$OMEGA 0.177 ;    18 BGIP

$SIGMA 0.0185 ; Paracetamol
$SIGMA 0.0123 ;    Glucose
$SIGMA 0.129 ;      GLP-1
$SIGMA 0.314 ;        GIP

$ESTIMATION METHOD=0 SIGL=6 NSIG=2 INTER MAXEVAL=9999 PRINT=1 NOABORT

;$ESTIMATION METHOD=IMP EONLY=1 ISAMPLE=1000 NITER=5
;$COVARIANCE UNCONDITIONAL MATRIX=R
;$SIM (43658) NSUB=1

;$TABLE ID TIME STUDY EVID CMT AMT RATE INSU BASI BW DV T2DM NOPRINT NOAPPEND ONEHEAD FILE=simdata126h
