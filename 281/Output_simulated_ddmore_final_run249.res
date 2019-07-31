05/02/2018 
11:43
$PROB B.dat 4-cRUN249
$INPUT ID TIME AMT=DOSE DV EVID MDV CMT RATE DLVL ROUT=DROP FORM=DROP
PERD=DROP SBP=DROP DBP ATEN=DROP METP=DROP BISO S1A2 D1A2=DROP
H1A2=DROP S3A4=DROP D3A4=DROP H3A4=DROP UPH=DROP GGT=DROP CK=DROP LDH
SGPT SGOT=DROP BIL SCR=DROP HR=DROP CARD=DROP SMK AGE=DROP HGHT=DROP
WGHT=DROP GEND=DROP BMI BSA=DROP LBM=DROP FAT=DROP CRCL


$DATA ../data/Simulated_Lid_B04_ddmore.csv IGNORE='C'
IGNORE=(TIME.LT.0)

$SUBR ADVAN5 TRANS1
$MODEL
COMP=(CENTRAL DEFDOSE DEFOBSERVATION)
COMP=(MEGX)
COMP=(GX)
COMP=(26XYL)

$PK
IF(AMT.GT.0)TOD=TIME
TAD=TIME-TOD

P1=1
P2=1
P3=1
P4=1
P5=1
P6=1
P7=1
IF(DLVL.GT.2)P1=0
IF(BIL.GT.0.53)P2=0
IF(LDH.GT.195)P3=0
IF(CRCL.LE.52.7)P4=0
IF(S1A2.EQ.3)P5=0
IF(SGPT.GT.11)P6=0
IF(BMI.GT.27.93)P7=0

TK12=THETA(1)
K12=TK12
TK23=THETA(2)
K23=TK23
TK14=THETA(3)
K14=TK14

T1K30=P1*THETA(4)+(1-P1)*THETA(5)
T2K30=T1K30+(1-P2)*THETA(6)
IF(T2K30.LE.0)T2K30=0.0001

T3K30=T2K30+(1-P4)*THETA(7)
IF(T3K30.LE.0)T3K30=0.0001

T4K30=T3K30+(1-P5)*THETA(8)
IF(T4K30.LE.0)T4K30=0.0001

T5K30=T4K30+(1-P7)*THETA(9)
IF(T5K30.LE.0)T5K30=0.0001

TK30=T5K30+(1-P6)*THETA(10)
IF(TK30.LE.0)TK30=0.0001


K30=TK30*EXP(ETA(1))

T1K40=P3*THETA(11)+(1-P3)*THETA(12)
TK40=T1K40+(1-P6)*THETA(13)
IF(TK40.LE.0)TK40=0.0001

K40=TK40*EXP(ETA(2))

TV1=P1*THETA(14)+(1-P1)*THETA(15)
V1=TV1*EXP(ETA(3))
TVM=THETA(16)
V2=TVM
V3=TVM
V4=TVM


S1=V1
S2=V2
S3=V3
S4=V4

$ERROR

IPRED=F
IRES=DV-IPRED

Q1=1
Q2=0
Q3=0

Y1=F+EPS(1)
Y2=F+EPS(2)
Y3=F+EPS(3)
Y4=F+EPS(4)

IF(CMT.EQ.2) THEN
  Q1=0
  Q2=1
  Q3=0
ENDIF

IF(CMT.EQ.3) THEN
  Q1=0
  Q2=0
  Q3=1
ENDIF

IF(CMT.EQ.4) THEN
  Q1=0
  Q2=0
  Q3=0
ENDIF

Y=Y1*Q1+Y2*Q2+Y3*Q3+Y4*(1-Q1)*(1-Q2)*(1-Q3)

W     = SQRT(SIGMA(1,1))       
IWRES = IRES/W			; IWRES for LID

$THETA
(0.03 FIX)     ;K12
(0,1)          ;K23
(0.007 FIX)    ;K14
(0,1)          ;K30  FOR DLVL <= 2
(0,1)          ;K30  FOR DLVL >  2

(-2,0.1,2)     ;MODIFICATION OF K30 FOR BIL > 0.53
(-2,0.1,2)     ;MODIFICATION OF K30 FOR CRCL <= 52.7
(-2,0.1,2)     ;MODIFICATION OF K30 FOR S1A2 = 3 (PRESENT)
(-2,0.1,2)     ;MODIFICATION OF K30 FOR BMI > 27.93
(-2,0.1,2)     ;MODIFICATION OF K30 FOR SGPT > 11

(0,1)          ;K40  FOR LDH <= 195
(0,1)          ;K40  FOR LDH >  195
(-2,0.1,2)     ;MODIFICATION OF K40 FOR SGPT > 11

(1,1000)       ;V1   FOR DLVL <= 2
(1,1000)       ;V1   FOR DLVL >  2
(100 FIX)      ;V2,3,4

$OMEGA
0.1            ;ETA ON K30
0.1            ;ETA ON K40
0.1            ;ETA ON V1

$SIGMA

10            ;EPS ADD. FOR LID 
10            ;EPS ADD. FOR MEGX
10            ;EPS ADD. FOR GX
10            ;EPS ADD. FOR 26-XYL  

$EST MAXEVAL=3000 PRINT=5 SIGDIG=3 NOABORT POSTHOC
$COV 

$TABLE MDV EVID ID TIME TAD DV CMT AMT RATE DLVL BIL CRCL S1A2 BMI SGPT LDH IPRED 
	IRES CWRES IWRES NPDE ETA1 ETA2 ETA3 NOPRINT FILE=ddmore_final_run249.tab ONEHEADER


  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
             
 (WARNING  79) SIGMA IS USED ON THE RIGHT. WITH A SUBSEQUENT RUN, IF AN
 INITIAL ESTIMATE OF A DIAGONAL BLOCK OF SIGMA IS TO BE COMPUTED BY
 NONMEM, THAT BLOCK WILL BE SET TO AN IDENTITY MATRIX DURING THAT
 COMPUTATION. THIS COULD LEAD TO AN ARITHMETIC EXCEPTION.*

 * THE MAXIMUM NUMBER OF WARNINGS OF ONE OR MORE TYPES WAS REACHED.
 IT IS POSSIBLE THAT SOME WARNING MESSAGES WERE SUPPRESSED.
  
License Registered to: BAST Inc. Ltd
Expiration Date:    14 JAN 2019
Current Date:        5 FEB 2018
Days until program expires : 344
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.3.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 B.dat 4-cRUN249
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:    16857
 NO. OF DATA ITEMS IN DATA SET:  18
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   4
 MDV DATA ITEM IS DATA ITEM NO.:  6
0INDICES PASSED TO SUBROUTINE PRED:
   5   2   3   8   0   0   7   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DOSE DV EVID MDV CMT RATE DLVL DBP BISO S1A2 LDH SGPT BIL SMK BMI CRCL
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 TAD IPRED IRES IWRES
0FORMAT FOR DATA:
 (E3.0,E18.0,E8.0,E19.0,3E1.0,E10.0,E4.0,E2.0,E6.0,E1.0/6E11.0)

 TOT. NO. OF OBS RECS:     1989
 TOT. NO. OF INDIVIDUALS:    325
0LENGTH OF THETA:  16
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   4
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.3000E-01     0.3000E-01     0.3000E-01
  0.0000E+00     0.1000E+01     0.1000E+07
  0.7000E-02     0.7000E-02     0.7000E-02
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
 -0.2000E+01     0.1000E+00     0.2000E+01
 -0.2000E+01     0.1000E+00     0.2000E+01
 -0.2000E+01     0.1000E+00     0.2000E+01
 -0.2000E+01     0.1000E+00     0.2000E+01
 -0.2000E+01     0.1000E+00     0.2000E+01
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
 -0.2000E+01     0.1000E+00     0.2000E+01
  0.1000E+01     0.1000E+04     0.1000E+07
  0.1000E+01     0.1000E+04     0.1000E+07
  0.1000E+03     0.1000E+03     0.1000E+03
0INITIAL ESTIMATE OF OMEGA:
 0.1000E+00
 0.0000E+00   0.1000E+00
 0.0000E+00   0.0000E+00   0.1000E+00
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+02
 0.0000E+00   0.1000E+02
 0.0000E+00   0.0000E+00   0.1000E+02
 0.0000E+00   0.0000E+00   0.0000E+00   0.1000E+02
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 RELATIVE TOLERANCE (TOL):                  -1
 ABSOLUTE TOLERANCE-ADVAN 9,13 ONLY (ATOL): -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 RESUME COV ANALYSIS (RESUME):               NO
0TABLES STEP OMITTED:    NO
 NO. OF TABLES:           1
 SEED NUMBER (SEED):    11456
 RANMETHOD:
 MC SAMPLES (ESEED):    300
 WRES SQUARE ROOT TYPE:            EIGENVALUE
0-- TABLE   1 --
 PRINTED:                NO
 HEADERS:               ONE
 FILE TO BE FORWARDED:   NO
 FORMAT:                S1PE11.4
 LFORMAT:
 RFORMAT:
0USER-CHOSEN ITEMS:
 MDV EVID ID TIME TAD DV CMT DOSE RATE DLVL BIL CRCL S1A2 BMI SGPT LDH IPRED IRES CWRES IWRES NPDE ETA1 ETA2 ETA3
1DOUBLE PRECISION PREDPP VERSION 7.3.0

 GENERAL LINEAR KINETICS MODEL (ADVAN5)
0MODEL SUBROUTINE USER-SUPPLIED - ID NO. 9999
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0RATE CONSTANT PARAMETERS - ASSIGNMENT OF ROWS IN GG
            TO COMPT.
  FROM      1    2    3    4    5
  COMPT.
    1       *    1    -    3    -
    2       -    *    2    -    -
    3       -    -    *    -    4
    4       -    -    -    *    5
             * LINK FROM A COMPARTMENT TO ITSELF IS NOT POSSIBLE
             - LINK BETWEEN THESE COMPARTMENTS IS NOT DEFINED FOR THIS MODEL
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         CENTRAL      ON         YES        YES        YES        YES
    2         MEGX         ON         YES        YES        NO         NO
    3         GX           ON         YES        YES        NO         NO
    4         26XYL        ON         YES        YES        NO         NO
    5         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            6           *           *           *           *
    2            7           *           *           *           *
    3            8           *           *           *           *
    4            9           *           *           *           *
    5            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      5
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   3
   DOSE RATE DATA ITEM IS DATA ITEM NO.:     8
   COMPT. NO. DATA ITEM IS DATA ITEM NO.:    7

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order

 ESTIMATION STEP OMITTED:                 NO  
 ANALYSIS TYPE:                           POPULATION
 EPS-ETA INTERACTION:                     NO  
 POP. ETAS OBTAINED POST HOC:             YES 
 NO. OF FUNCT. EVALS. ALLOWED:            3000
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES 
 ESTIMATE OUTPUT TO MSF:                  NO  
 ABORT WITH PRED EXIT CODE 1:             NO  
 IND. OBJ. FUNC. VALUES SORTED:           NO  
 NUMERICAL DERIVATIVE 
       FILE REQUEST (NUMDER):             NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP): 0           
 ETA HESSIAN EVALUATION METHOD (ETADER):  0           
 INITIAL ETA FOR MAP ESTIMATION (MCETA):  0           
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):    100         
 GRADIENT SIGDIGITS OF 
       FIXED EFFECTS PARAMETERS (SIGL):   100         
 EXCLUDE TITLE (NOTITLE):                 NO 
 EXCLUDE COLUMN LABELS (NOLABEL):         NO 
 NOPRIOR SETTING (NOPRIOR):               OFF
 NOCOV SETTING (NOCOV):                   OFF
 DERCONT SETTING (DERCONT):               OFF
 ABSOLUTE TOLERANCE-ADVAN 9,13 ONLY(ATOL):-100        
 FINAL ETA RE-EVALUATION (FNLETA):        ON 
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS 
       IN SHRINKAGE (ETASTYPE):           NO 
 NON-INFL. ETA CORRECTION (NONINFETA):    OFF
 FORMAT FOR ADDITIONAL FILES (FORMAT):    S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):     TSOL
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:  NO
 EM OR BAYESIAN METHOD USED:                NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=NPRED
 RES=NRES
 WRES=NWRES
 IWRS=NIWRES
 IPRD=NIPRED
 IRS=NIRES

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:   26680.1140072524        NO. OF FUNC. EVALS.:  19
 CUMULATIVE NO. OF FUNC. EVALS.:       19
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E+00  1.0000E+00
             1.0000E-01  1.0000E+03  1.0000E+03  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E+01  1.0000E+01  1.0000E+01  1.0000E+01
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:  -8.1663E+03  9.7147E+02 -7.5886E+02  0.0000E+00  0.0000E+00 -1.1952E+02  0.0000E+00  0.0000E+00  5.6794E+02  0.0000E+00
             0.0000E+00  3.7533E+03 -2.5732E+03 -8.3485E+02 -1.4181E+02 -2.0081E+03 -2.3569E+04 -8.5924E+03 -3.0585E+03 -9.3943E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:   11997.7578471038        NO. OF FUNC. EVALS.:  20
 CUMULATIVE NO. OF FUNC. EVALS.:      120
 NPARAMETR:  1.2307E+00  9.8974E-01  1.2617E+00  1.0000E-01  1.0000E-01  2.7568E-01  1.0000E-01  1.0000E-01  5.8586E-01  1.0000E+00
             1.0000E-01  1.0844E+03  1.2533E+03  2.0153E-01  1.3036E-01  2.0993E-01  7.1538E+01  2.8405E+01  2.5471E+01  7.3150E+00
 PARAMETER:  3.0758E-01  8.9685E-02  3.3249E-01  1.0000E-01  1.0000E-01  2.7736E-01  1.0000E-01  1.0000E-01 -4.3468E-01  1.0000E-01
             1.0000E-01  1.8108E-01  3.2602E-01  4.5038E-01  2.3255E-01  4.7080E-01  1.0838E+00  6.2200E-01  5.6748E-01 -5.6329E-02
 GRADIENT:  -1.3999E+03  2.6015E+02 -2.1107E+02  0.0000E+00  0.0000E+00 -4.2726E+01  0.0000E+00  0.0000E+00 -4.0656E+00  0.0000E+00
             0.0000E+00  5.1576E+02 -3.0384E+02 -2.3305E+02 -3.8341E+01 -4.0681E+02 -2.8961E+03 -1.5475E+03 -7.0815E+02  4.6951E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:   11954.0065109054        NO. OF FUNC. EVALS.:  21
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  1.2951E+00  1.0713E+00  1.4537E+00  1.0000E-01  1.0000E-01  4.4500E-01  1.0000E-01  1.0000E-01  6.6416E-01  1.0000E+00
             1.0000E-01  1.0593E+03  1.2581E+03  2.9459E-01  2.7909E-01  2.0489E-01  7.0081E+01  2.7060E+01  2.4786E+01  5.9534E+00
 PARAMETER:  3.5859E-01  1.6887E-01  4.7409E-01  1.0000E-01  1.0000E-01  4.5248E-01  1.0000E-01  1.0000E-01 -3.0924E-01  1.0000E-01
             1.0000E-01  1.5765E-01  3.2978E-01  6.4021E-01  6.1318E-01  4.5865E-01  1.0735E+00  5.9774E-01  5.5384E-01 -1.5931E-01
 GRADIENT:  -1.1945E+03  2.7593E+02 -1.5980E+02  0.0000E+00  0.0000E+00 -4.0684E+00  0.0000E+00  0.0000E+00  8.3800E+01  0.0000E+00
             0.0000E+00  4.7816E+02 -2.8536E+02 -1.7020E+02 -2.1894E-01 -4.1478E+02 -2.9452E+03 -1.5360E+03 -6.8892E+02 -6.3522E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:   10882.6415379516        NO. OF FUNC. EVALS.:  20
 CUMULATIVE NO. OF FUNC. EVALS.:      328
 NPARAMETR:  1.4346E+00  9.0777E-01  1.0211E+00  1.0000E-01  1.0000E-01  8.0361E-01  1.0000E-01  1.0000E-01  6.4541E-01  1.0000E+00
             1.0000E-01  1.2677E+03  1.8692E+03  4.0655E-01  2.8140E-01  1.7761E-01  1.5195E+02  4.1270E+01  3.5630E+01  6.4005E+00
 PARAMETER:  4.6090E-01  3.2393E-03  1.2092E-01  1.0000E-01  1.0000E-01  8.5151E-01  1.0000E-01  1.0000E-01 -3.3787E-01  1.0000E-01
             1.0000E-01  3.3742E-01  7.2597E-01  8.0127E-01  6.1730E-01  3.8720E-01  1.4605E+00  8.0877E-01  7.3530E-01 -1.2311E-01
 GRADIENT:  -4.8569E+02  1.1960E+02 -2.1593E+02  0.0000E+00  0.0000E+00 -6.2870E+00  0.0000E+00  0.0000E+00  4.0463E+01  0.0000E+00
             0.0000E+00  2.8322E+02  1.0849E+02 -1.0725E+01  1.6853E+01 -2.8638E+02 -1.2555E+03 -5.5195E+02 -3.0757E+02 -5.9695E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:   10399.4391460508        NO. OF FUNC. EVALS.:  20
 CUMULATIVE NO. OF FUNC. EVALS.:      433
 NPARAMETR:  1.8169E+00  1.0140E+00  2.0651E+00  1.0000E-01  1.0000E-01  8.2084E-01  1.0000E-01  1.0000E-01  5.9937E-01  1.0000E+00
             1.0000E-01  1.1984E+03  1.7369E+03  6.0552E-01  2.0516E-01  2.9670E-01  4.0881E+02  6.0255E+01  4.8721E+01  6.4912E+00
 PARAMETER:  6.9715E-01  1.1386E-01  8.2517E-01  1.0000E-01  1.0000E-01  8.7215E-01  1.0000E-01  1.0000E-01 -4.1188E-01  1.0000E-01
             1.0000E-01  2.8118E-01  6.5251E-01  1.0005E+00  4.5931E-01  6.4378E-01  1.9553E+00  9.9800E-01  8.9176E-01 -1.1607E-01
 GRADIENT:  -2.9549E+01 -8.3377E+00  8.9724E+00  0.0000E+00  0.0000E+00  1.1930E+00  0.0000E+00  0.0000E+00 -4.6727E+00  0.0000E+00
             0.0000E+00  6.2047E+00 -7.0033E-01  1.1247E+01 -1.8982E+00  2.6953E+00  1.8110E+00 -1.2812E+00  2.7712E+01 -8.1745E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:   10398.0071906975        NO. OF FUNC. EVALS.:  20
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  1.8787E+00  9.9192E-01  1.9492E+00  1.0000E-01  1.0000E-01  7.5826E-01  1.0000E-01  1.0000E-01  6.0320E-01  1.0000E+00
             1.0000E-01  1.1816E+03  1.7258E+03  5.3262E-01  2.1078E-01  2.8395E-01  4.0922E+02  6.0273E+01  4.6857E+01  6.4884E+00
 PARAMETER:  7.3057E-01  9.1892E-02  7.6743E-01  1.0000E-01  1.0000E-01  7.9800E-01  1.0000E-01  1.0000E-01 -4.0550E-01  1.0000E-01
             1.0000E-01  2.6705E-01  6.4609E-01  9.3632E-01  4.7281E-01  6.2182E-01  1.9558E+00  9.9815E-01  8.7226E-01 -1.1628E-01
 GRADIENT:  -7.8334E-01  2.0025E-01 -1.8507E-01  0.0000E+00  0.0000E+00 -1.4986E-01  0.0000E+00  0.0000E+00  1.1338E-01  0.0000E+00
             0.0000E+00 -6.1249E-01  3.1298E-01 -1.4792E-01 -1.2735E-02  1.0770E-01  1.1126E+00 -4.6959E-03 -3.4019E-01 -2.7296E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:   10398.0047151530        NO. OF FUNC. EVALS.:  37
 CUMULATIVE NO. OF FUNC. EVALS.:      668
 NPARAMETR:  1.8812E+00  9.9202E-01  1.9516E+00  1.0000E-01  1.0000E-01  7.6318E-01  1.0000E-01  1.0000E-01  6.0330E-01  1.0000E+00
             1.0000E-01  1.1828E+03  1.7255E+03  5.3376E-01  2.1083E-01  2.8388E-01  4.0948E+02  6.0335E+01  4.6921E+01  6.4894E+00
 PARAMETER:  7.3193E-01  9.1984E-02  7.6865E-01  1.0000E-01  1.0000E-01  8.0376E-01  1.0000E-01  1.0000E-01 -4.0534E-01  1.0000E-01
             1.0000E-01  2.6806E-01  6.4593E-01  9.3739E-01  4.7295E-01  6.2168E-01  1.9562E+00  9.9866E-01  8.7294E-01 -1.1621E-01
 GRADIENT:   3.2596E-02 -5.1250E-02  2.6761E-02  0.0000E+00  0.0000E+00 -2.4326E-03  0.0000E+00  0.0000E+00  1.5800E-02  0.0000E+00
             0.0000E+00 -3.3339E-02  1.0834E-02  1.4852E-02 -4.9993E-03  5.1652E-02 -4.0174E-02  2.8752E-02 -5.9639E-02 -1.3201E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:   10398.0047151530        NO. OF FUNC. EVALS.:   0
 CUMULATIVE NO. OF FUNC. EVALS.:      668
 NPARAMETR:  1.8812E+00  9.9202E-01  1.9516E+00  1.0000E-01  1.0000E-01  7.6318E-01  1.0000E-01  1.0000E-01  6.0330E-01  1.0000E+00
             1.0000E-01  1.1828E+03  1.7255E+03  5.3376E-01  2.1083E-01  2.8388E-01  4.0948E+02  6.0335E+01  4.6921E+01  6.4894E+00
 PARAMETER:  7.3193E-01  9.1984E-02  7.6865E-01  1.0000E-01  1.0000E-01  8.0376E-01  1.0000E-01  1.0000E-01 -4.0534E-01  1.0000E-01
             1.0000E-01  2.6806E-01  6.4593E-01  9.3739E-01  4.7295E-01  6.2168E-01  1.9562E+00  9.9866E-01  8.7294E-01 -1.1621E-01
 GRADIENT:   3.2596E-02 -5.1250E-02  2.6761E-02  0.0000E+00  0.0000E+00 -2.4326E-03  0.0000E+00  0.0000E+00  1.5800E-02  0.0000E+00
             0.0000E+00 -3.3339E-02  1.0834E-02  1.4852E-02 -4.9993E-03  5.1652E-02 -4.0174E-02  2.8752E-02 -5.9639E-02 -1.3201E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      668
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
 #TERE:
 Elapsed estimation time in seconds:   583.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance time in seconds:   589.18
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    10398.005       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12     
          TH13      TH14      TH15      TH16     
 
         3.00E-02  1.88E+00  7.00E-03  9.92E-01  1.95E+00  1.00E-01  1.00E-01  7.63E-01  1.00E-01  1.00E-01  6.03E-01  1.00E+00
          1.00E-01  1.18E+03  1.73E+03  1.00E+02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


            ETA1      ETA2      ETA3   
 
 ETA1
+        5.34E-01
 
 ETA2
+        0.00E+00  2.11E-01
 
 ETA3
+        0.00E+00  0.00E+00  2.84E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


            EPS1      EPS2      EPS3      EPS4   
 
 EPS1
+        4.09E+02
 
 EPS2
+        0.00E+00  6.03E+01
 
 EPS3
+        0.00E+00  0.00E+00  4.69E+01
 
 EPS4
+        0.00E+00  0.00E+00  0.00E+00  6.49E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


            ETA1      ETA2      ETA3   
 
 ETA1
+        7.31E-01
 
 ETA2
+        0.00E+00  4.59E-01
 
 ETA3
+        0.00E+00  0.00E+00  5.33E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


            EPS1      EPS2      EPS3      EPS4   
 
 EPS1
+        2.02E+01
 
 EPS2
+        0.00E+00  7.77E+00
 
 EPS3
+        0.00E+00  0.00E+00  6.85E+00
 
 EPS4
+        0.00E+00  0.00E+00  0.00E+00  2.55E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                                   FIRST ORDER                                  ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12  
             TH13      TH14      TH15      TH16      OM11      OM12      OM13      OM22      OM23      OM33      SG11      SG12  
            SG13      SG14      SG22      SG23      SG24      SG33      SG34      SG44  
 
 TH 1
+       .........
 
 TH 2
+       .........  1.15E+02
 
 TH 3
+       ......... ......... .........
 
 TH 4
+       .........  2.32E-01 .........  3.21E+02
 
 TH 5
+       .........  6.80E-02 ......... -3.41E+01  4.02E+01
 
 TH 6
+       .........  0.00E+00 .........  0.00E+00  0.00E+00  0.00E+00
 
 TH 7
+       .........  0.00E+00 .........  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       ......... -3.62E-02 .........  1.82E+01  3.69E+00  0.00E+00  0.00E+00  2.19E+01
 
 TH 9
+       .........  0.00E+00 .........  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       .........  0.00E+00 .........  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       .........  0.00E+00 .........  9.37E-06  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  1.70E+03
 
 TH12
+       .........  0.00E+00 .........  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH13
+       .........  0.00E+00 .........  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
          0.00E+00
 
 TH14
+       .........  5.21E-10 ......... -7.23E-09 -4.78E-10  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00 -2.93E-09  0.00E+00
          0.00E+00  2.48E-04
 
 TH15
+       .........  0.00E+00 .........  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
          0.00E+00 -2.50E-05  7.23E-05
 
 TH16
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM11
+       .........  6.74E-02 ......... -1.18E+02 -2.02E+01  0.00E+00  0.00E+00 -1.39E+01  0.00E+00  0.00E+00  1.68E-11  0.00E+00
          0.00E+00  0.00E+00  0.00E+00 .........  1.52E+02
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      TH12  
             TH13      TH14      TH15      TH16      OM11      OM12      OM13      OM22      OM23      OM33      SG11      SG12  
            SG13      SG14      SG22      SG23      SG24      SG33      SG34      SG44  
 
 OM22
+       ......... -4.64E-12 .........  9.15E-11 -7.01E-12  0.00E+00  0.00E+00  6.86E-12  0.00E+00  0.00E+00 -4.75E+02  0.00E+00
          0.00E+00  0.00E+00  0.00E+00 .........  0.00E+00 ......... .........  7.15E+02
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       .........  0.00E+00 .........  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
          0.00E+00 -7.72E-02 -5.36E-02 .........  0.00E+00 ......... .........  0.00E+00 .........  3.44E+02
 
 SG11
+       .........  0.00E+00 .........  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  5.80E-10  0.00E+00
          0.00E+00 -9.32E-05 -8.03E-06 .........  0.00E+00 ......... .........  0.00E+00 .........  2.04E-01  1.08E-03
 
 SG12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 SG13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 SG14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG22
+       ......... -1.26E-03 .........  1.90E-08  1.26E-09  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  7.72E-09  0.00E+00
          0.00E+00 -5.96E-12 -1.69E-12 .........  1.89E-09 ......... ......... -9.46E-09 ......... -5.35E-09  0.00E+00 .........
        ......... .........  6.51E-02
 
 SG23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... ......... .........
 
 SG24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... ......... ......... .........
 
 SG33
+       ......... -2.46E-03 ......... -2.44E-01 -5.89E-02  0.00E+00  0.00E+00 -1.47E-01  0.00E+00  0.00E+00  0.00E+00  0.00E+00
          0.00E+00  0.00E+00  0.00E+00 .........  4.64E-01 ......... .........  0.00E+00 .........  0.00E+00  0.00E+00 .........
        ......... .........  0.00E+00 ......... .........  7.87E-02
 
 SG34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... ......... ......... ......... ......... .........
 
 SG44
+       .........  0.00E+00 ......... -1.52E-06  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00 -1.22E+01  0.00E+00
          0.00E+00  4.76E-10  0.00E+00 .........  0.00E+00 ......... .........  1.56E+01 .........  4.27E-07  9.41E-11 .........
        ......... ......... -1.25E-09 ......... .........  1.84E-09 .........  4.43E+00
 
 #CPUT: Total CPU Time in Seconds,     1505.394
Stop Time: 
05/02/2018 
12:08
