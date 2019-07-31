

$PROB Logit model with KPD approach EMAX model covariate CRCL
$INPUT C	ID	WEEK=TIME	MDV	EVID	GHFS=DV	AMT	ADDL	II	AGE	SEX	WT	HT	SERUM SERCORR	WTCOR	SEXCOEFF	BSA	CLCR
$DATA ../Datasets/Simulated_HFSmodel.csv IGNORE=C
$SUBROUTINES ADVAN1 TRANS1 

$PK
TCLCR=THETA(12)
TVK=THETA(1)

K=TVK*EXP(ETA(1))


IF (TIME.EQ.0) BCLCR=CLCR
IIV=ETA(2)+(BCLCR-75.5)*TCLCR

B00=THETA(2)                    
B10=THETA(3)                    
B20=THETA(4)                    
EMAX0=THETA(5)
ED50=THETA(6)  

B01=THETA(7)
B11=THETA(8)
B21=THETA(9)

EMAX1=THETA(10)
EMAX2=THETA(11)

IF (TIME.EQ.0) SWM1=0


$ERROR
CALLFL=0
IPRED=F
EMAX=EMAX0 
IF(SWM1.EQ.1) EMAX=EMAX1
IF(SWM1.EQ.2) EMAX=EMAX2
EFF=EMAX*(F*K)/((F*K)+(ED50))
A0=0
A1=0
SPREC=SWM1
IF (SWM1.EQ.0) THEN
        A0=B00
        A1=A0+B01
ENDIF
IF (SWM1.EQ.1) THEN
        A0=B10
        A1=A0+B11
ENDIF
IF (SWM1.EQ.2) THEN
        A0=B20
        A1=A0+B21
ENDIF
A0=A0-EFF+IIV
A1=A1-EFF+IIV

SWM1=GHFS

PC0=EXP(A0)/(1+EXP(A0))
PC1=EXP(A1)/(1+EXP(A1))
PC2=1


P0=PC0
P1=PC1-PC0
P2=PC2-PC1


Y=-1
   IF (DV.LT.0.5)  Y=P0
   IF (DV.GE.0.5.AND.DV.LT.1.5)  Y=P1
   IF (DV.GE.1.5.AND.DV.LT.2.5)  Y=P2

$THETA
      (0,0.159)                  ; 1 TVK
      (4.62)                     ; 2 B00
      (0.683)                    ; 3 B10
      (1.99)                     ; 4 B20
      (3.8)                      ; 5 EMAX0
      (0,13000)                  ; 6 ED50
      (0,0.602)                  ; 7* B01
      (0,5.24)                   ; 8* B11
      (0,0.322)                  ; 9 B21
      (0,6.3)                    ;10 EMAX1
      (0,9.9)                    ;11 EMAX2
      (0,0.00552)                ;12 TCLCR


$OMEGA BLOCK(2)
       0.468
       0.402 0.8 

$EST METHOD=1 MAXEVALS=9999 PRINT=5 LIKE LAPLACE SIGDIGITS=1 SLOW
NOABORT 
$TABLE ID TIME AMT MDV SPREC IPRED K TVK TCLCR
P0 P1 P2 A0 A1 PC0 PC1 PC2 ETA1 ETA2 EMAX ED50 EFF
ONEHEADER NOPRINT FILE=ddmore159.tab
$COV SLOW UNCONDITIONAL MATRIX=S

