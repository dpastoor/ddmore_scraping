; The Multistate Tuberculosis Pharmacometric (MTP) model combined with 
; the General Pharmacodynamic Interaction model (GPDI)
; Isoniazid adaptive resistance
; M. tuberculosis B1585 (cfu/ml)
; Data with rifampicin (RIF)(R), isoniazid (INH)(H) and 
; ethambutol (EMB)(E) in mono and combination exposure (mg/L)
; All MTP model parameters were fixed in this model when evaluating
; the PD interactions
; NM version 7.3
$PROBLEM    MTP-GPDI RIF INH EMB MTB1585
$INPUT      ID REPL TIME EVID DV RIF INH EMB PZA SM PAS EXPR PH THRP
            LOQ ARG NDV EVIDLQ NATG IDEXPR PLOT
$DATA      ../Data/Erasmus_MC_Beijing-1585_Drug_data_and_NATG-2016-10-05.csv
            IGNORE=@ ACCEPT=(NATG.EQ.1) ACCEPT=(EXPR.EQ.1)
            ACCEPT=(EXPR.EQ.2) ACCEPT=(EXPR.EQ.3) ACCEPT=(EXPR.EQ.7)
            ACCEPT=(EXPR.EQ.8) ACCEPT=(EXPR.EQ.11) ACCEPT=(EXPR.EQ.13)
$SUBROUTINE ADVAN13 TOL=9
$MODEL      NCOMP=8 COMP=(FBUGS) COMP=(SBUGS) COMP=(NBUGS) COMP=(INH)
            COMP=(ARON) COMP=(AROFF) COMP=(RIF) COMP=(EMB)
$PK
; MTP model                                                                
TVKG=THETA(1)           ; Growth rate, F
KFSLIN=THETA(2)/100     ; Rate parameter, KFS, time dependent
KFN=THETA(3)/1000000    ; Transfer rate F -> N
KSF=THETA(4)/10         ; Transfer rate S -> F
KSN=THETA(5)            ; Transfer rate S -> N
KNS=THETA(6)/100        ; Transfer rate N -> S

TVF0=THETA(7)*1000      ; Initial F bacterial number with IIV
TVS0=THETA(8)*1000      ; Initial S bacterial number

; INH(H) exposure-response relationship (estimated with mono data)
HFDEMAX=THETA(9)        ; EMAX FD
HFDEC50=THETA(10)       ; EC50 FD
HFDGAM=THETA(11)        ; GAMMA FD
HSDEMAX=THETA(12)       ; EMAX SD
HSDEC50=THETA(13)       ; EC50 SD
HSDGAM=THETA(14)        ; GAMMA for SD

; INH(H) adaptive resistance (estimated with mono data)
KON=THETA(15)           ; kON
KOFF=THETA(16)          ; kOFF
ARLINFD=THETA(17)       ; ARLIN FD (Slope)
ARLINSD=THETA(18)       ; ARLIN SD (Slope)

; RIF(R) exposure-response relationship (estimated with mono data)
RFDEMAX=THETA(19)       ; EMAX on FD
RFDEC50=THETA(20)       ; EC50 on FD
RFGEMAX=THETA(21)       ; EMAX on FG
RFGEC50=THETA(22)       ; EC50 on FG
RFGGAM=THETA(23)        ; GAMMA on FG
RSDEMAX=THETA(24)       ; EMAX for SD
RSDEC50=THETA(25)       ; EC50 for SD
RNDK=THETA(26)          ; k for ND

; EMB(E) exposure-response relationship (estimated with mono data)
EFDEMAX=THETA(27)       ; EMAX on FD
EFDEC50=THETA(28)       ; EC50 on FD
EFDGAM=THETA(29)        ; GAMMA on FD
ESDK=THETA(30)          ; k on SD

; GPDI model (estimated with all data and MTP model fixed)
; INH(H) - RIF(R) PD interaction
FDIRH=THETA(31)         ; FD interaction RIF-INH
FDIHR=THETA(32)         ; FD interaction INH-RIF
SDIRH=THETA(33)         ; SD interaction RIF-INH
SDIHR=THETA(34)         ; SD interaction INH-RIF

; INH(H) - EMB(E) PD interaction
FDIEH=THETA(35)         ; FD interaction EMB-INH
FDIHE=THETA(36)         ; FD interaction INH-EMB
IF(A(8).GT.0) THEN
SDIEH=THETA(37)         ; SD interaction EMB-INH
ELSE
SDIEH=0
ENDIF
IF(A(4).GT.0) THEN
SDIHE=THETA(38)         ; SD interaction INH-EMB
ELSE
SDIHE=0                 ; SD interaction INH-EMB
ENDIF

; RIF(R) - EMB(E) PD interaction
IF(A(8).GT.0) THEN
FDIER=THETA(39)         ; FD interaction EMB-RIF
ELSE
FDIER=0
ENDIF
FDIRE=THETA(40)         ; FD interaction RIF-EMB
IF(A(8).GT.0) THEN
SDIER=THETA(41)         ; SD interaction EMB-RIF
ELSE
SDIER=0
ENDIF
SDIRE=THETA(42)         ; SD interaction RIF-EMB

; EMB(E) - (RIF(R) - INH(H)) PD interaction
IF(A(8).GT.0) THEN
SDIERH=THETA(43)        ; SD interaction EMB-(RIF-INH)
ELSE
SDIERH=0
ENDIF

KG=TVKG
F0=TVF0
S0=TVS0

A_0(1)=F0                         ; Initial F bacterial number with IIV
A_0(2)=S0                         ; Initial S bacterial number
A_0(3)=0.00001                    ; Initial N bacterial number
A_0(4)=INH                        ; INH concentration
A_0(5)=0                          ; ARON
A_0(6)=1                          ; AROFF
A_0(7)=RIF                        ; RIF concentration
A_0(8)=EMB                        ; EMB concentration

$DES
; MTP model                                                                       
GROWTHFUNC=KG                     ; Exponential growth function
IF(GROWTHFUNC.LT.0) GROWTHFUNC=0  ; Keep GROWTHFUNC from turning negative

KFS=KFSLIN*T                      ; Time dependent linear transfer, F -> S

; INH adaptive resistance
AREFD=1+(ARLINFD*A(5))            ; Linear adaptive resistance function FD
ARESD=1+(ARLINSD*A(5))            ; Linear adaptive resistance function SD

HFDEC50A=HFDEC50*AREFD
HSDEC50A=HSDEC50*ARESD

; INH (H) drug effect on F
INHFD=1*A(4)**HFDGAM/((HFDEC50A*(1+(FDIRH*A(7)/(RFDEC50+A(7))))*(1+(FDIEH*A(8)/(EFDEC50+A(8)))))**HFDGAM+A(4)**HFDGAM)

; INH (H) drug effect on S
INHSD=HSDEMAX*A(4)**HSDGAM/((HSDEC50A*(1+(SDIRH*(1+(SDIERH))*A(7)/(RSDEC50+A(7))))*(1+(SDIEH)))**HSDGAM+A(4)**HSDGAM)

; RIF (R) drug effect on F
RIFFD=(RFDEMAX/HFDEMAX)*A(7)/((RFDEC50*(1+(FDIHR*A(4)/(HFDEC50A+A(4))))*(1+(FDIER)))+A(7))

; RIF (R) drug effect on FG
RIFFG=1-(RFGEMAX*A(7)**RFGGAM/(RFGEC50**RFGGAM+A(7)**RFGGAM))
IF(RIFEFG.LT.0) RIFEFG=0 ; Keep RIFEFG from turning negative

; RIF (R) drug effect on S
RIFSD=RSDEMAX*A(7)/((RSDEC50*(1+(SDIHR*A(4)/(HSDEC50A+A(4))))*(1+(SDIER)))+A(7))

; RIF (R) drug effect on N
RIFND=RNDK*A(7)

; EMB (E) drug effect on F
EMBFD=(EFDEMAX/HFDEMAX)*A(8)**EFDGAM/((EFDEC50*(1+(FDIHE*A(4)/(HFDEC50A+A(4))))*(1+(FDIRE*A(7)/(RFDEC50+A(7)))))**EFDGAM+A(8)**EFDGAM)

; EMB (E) drug effect on S
EMBSD=A(8)*(ESDK/((1+SDIHE*A(4))*(1+SDIRE*A(7)/RSDEC50+A(7))))

; Total drug effect with Bliss independence type interaction on F
FD=(INHFD+RIFFD+EMBFD-INHFD*RIFFD-INHFD*EMBFD-RIFFD*EMBFD+INHFD*RIFFD*EMBFD)*HFDEMAX

; Total drug effect on FG (only RIF was found to have effect)
FG=RIFFG

; Total drug effect with Bliss independence type interaction on S
SD=INHSD+RIFSD+EMBSD

; Total drug effect on N (only RIF was found to have effect)
ND=RIFND

DADT(1)=A(1)*FG*GROWTHFUNC+KSF*A(2)-KFS*A(1)-KFN*A(1)-FD*A(1)   ; F
DADT(2)=KFS*A(1)+KNS*A(3)-KSN*A(2)-KSF*A(2)-SD*A(2)             ; S
DADT(3)=KSN*A(2)+KFN*A(1)-KNS*A(3)-ND*A(3)                      ; N
DADT(4)=0                                                       ; INH
DADT(5)=KON*A(6)*A(4)-KOFF*A(5)                                 ; ARON
DADT(6)=KOFF*A(5)-KON*A(6)*A(4)                                 ; AROFF
DADT(7)=0                                                       ; RIF
DADT(8)=0                                                       ; EMB

$ERROR                                                                       
FBUGS=A(1)             ; F
SBUGS=A(2)             ; S
NBUGS=A(3)             ; N
TOTBUGS=A(1)+A(2)+A(3) ; F+S+N
CINH=A(4)              ; INH concentration
ARON=A(5)              ; Adaptive resistance on
AROFF=A(6)             ; Adaptive resistance off
CRIF=A(7)              ; RIF concentration
CEMB=A(8)              ; EMB concentration

F_FLAG=0
IPRED=LOG(A(1)+A(2))
IRES=DV-IPRED
STD=SQRT(SIGMA(1))     ; Standard deviation for additive residual error on log scale
IWRES=IRES/STD
Y=IPRED+EPS(1)

; LOQ
LLOQ=1.6094

; M3 method
IF (DV.LT.LLOQ) THEN
F_FLAG=1
MDVRES=1
DUM=(LLOQ-IPRED)/(SQRT(SIGMA(1)))
CUMD=PHI(DUM)
Y=CUMD
ENDIF

$THETA
; MTP model
(0,0.796) FIX     ; 1 kG
(0,0.166) FIX     ; 2 kFSLIN (/100)
(0,0.897) FIX     ; 3 kFN (/1000000)
(0,0.145) FIX     ; 4 kSF (/10)
(0,0.186) FIX     ; 5 kSN
(0,0.123) FIX     ; 6 kNS (/100)
(0,209) FIX       ; 7 F0 (*1000)
(0,324) FIX       ; 8 S0 (*1000)
; INH(H) exposure-response relationship (estimated with mono data)
(0,22.2209) FIX   ; 9 HFDEMAX
(0,0.168431) FIX  ; 10 HFDEC50
(0,1.90157) FIX   ; 11 HFDGAMMA
(0,8.55316) FIX   ; 12 HSDEMAX
(0,0.0328672) FIX ; 13 HSDEC50
(0,1.74098) FIX   ; 14 HSDGAMMA
; INH(H) adaptive resistance (estimated with mono data)
(0,0.0205994) FIX ; 15 kON
0 FIX             ; 16 kOFF
(0,522.42) FIX    ; 17 ARLIN FD
(0,2352.28) FIX   ; 18 ARLIN SD
; RIF(R) exposure-response relationship (estimated with mono data)
(0,1.96922) FIX   ; 19 RFDEMAX
(0,0.0030258) FIX ; 20 RFDEC50
(0,1) FIX         ; 21 RFGEMAX
(0,0.388407) FIX  ; 22 RFGEC50
(0,2.80234) FIX   ; 23 RFGGAM
(0,1.79211) FIX   ; 24 RSDEMAX
(0,0.0112798) FIX ; 25 RSDEC50
(0,3.28587) FIX   ; 26 RNDK
; EMB(E) exposure-response relationship (estimated with mono data)
(0,2.2073) FIX    ; 27 EFDEMAX
(0,0.860332) FIX  ; 28 EFDEC50
(0,2.45751) FIX   ; 29 EFDGAMMA
(0,4.38781) FIX   ; 30 ESDk
; GPDI model
; INH(H) - RIF(R) PD interaction 
(-1,-0.683351)    ; 31 FDIRH
(-1,0) FIX        ; 32 FDIHR
(-1,1.52908)      ; 33 SDIRH
(-1,10.7494)      ; 34 SDIHR
; INH(H) - EMB(E) PD interaction
(-1,1.80936)      ; 35 FDIEH
(-1,0) FIX        ; 36 FDIHE
(-1,0.0854746)    ; 37 SDIEH
(-1,91.4222)      ; 38 SDIHE
; RIF(R) - EMB(E) PD interaction
(-1,-0.662812)    ; 39 FDIER
(-1,-0.99999) FIX ; 40 FDIRE
(-1,1.70602)      ; 41 SDIER
(-1,479.458)      ; 42 SDIRE
; EMB(E) - (RIF(R) - INH(H)) PD interaction
(-1,-0.677036)    ; 43 SDIERH
$OMEGA 0 FIX ;
$SIGMA  0.936573  ; variance for add residual error on log scale
$ESTIMATION METHOD=1 LAPLACIAN INTER MAXEVAL=9999 NSIG=3 SIGL=9
            NOABORT
$COVARIANCE PRINT=E
$TABLE      ID TIME IPRED IRES IWRES CWRES NDV TOTBUGS FBUGS SBUGS
            NBUGS NATG LOQ CINH CRIF CEMB PLOT FD RIFFD INHFD
            EMBFD FG RIFFG SD RIFSD INHSD EMBSD ND RIFND EVID
            EXPR ONEHEADER NOPRINT FILE=sdtabMTPGPDI
$TABLE      ID TIME GROWTHFUNC KG KFN KFS KFSLIN KSF KSN KNS F0 S0 ETA(1) NATG CINH CRIF CEMB PLOT ONEHEADER
            NOPRINT FILE=patabMTPGPDI
$TABLE      ID TIME ONEHEADER NOPRINT FILE=cotabMTPGPDI
$TABLE      ID TIME ONEHEADER NOPRINT FILE=catabMTPGPDI
