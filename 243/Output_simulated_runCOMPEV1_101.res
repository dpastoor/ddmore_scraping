20/06/2017 
15:01


$PROBLEM Log_normal model. runCOMPEV1_101.    

$INPUT ID AGE NEUT PRE_TRE MAX_LEG AUC CMAX DBLOCK TIME DV DVID EVID
AGE_CAT NEUT_CAT PRE_TR_C MAX_CAT AUC_CAT CMAX_CAT

$DATA Simulated_event_data.csv
ACCEPT=(DVID.EQ.3,DVID.EQ.0)                                                        ; Select which DVID to accept

$SUBROUTINE ADVAN13 TRANS1 TOL=9

$MODEL NCOMPARTMENTS=1
       COMP=(CPT1)      	; TTE model			




$PK 

CENSORING = 2

IF (NEWIND.LE.1) SURVZ  = 1   		                ; Survival(0)=1 [we use this variable for storing survival function at start of observation interval for interval censored data]


;;;;;;;;;;;;;;;;;;;;; Start of TTE model ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

lambda      = THETA(1)*exp(ETA(1))
alpha	    = THETA(2)


pi  = 3.1415927
fac = (2*pi)**0.5

DEL = 1E-8


VAL          = EXP((THETA(3)/1000)*(AGE-55))


$DES 

num = log((T+DEL)/alpha)

pdf = (1.0/(fac*lambda*(T+DEL)))*exp(-(num**2)/(2*lambda*lambda))
DADT(1)  = VAL*pdf/(1-phi(num/lambda))               			            



$ERROR

CHAZ = A(1)				; cumulative hazard
SURV = EXP(-CHAZ)  			; probability of surviving to or beyond current time

num1 = log((TIME+DEL)/alpha)

pdf1 = (1.0/(fac*lambda*(TIME+DEL)))*exp(-(num1**2)/(2*lambda*lambda))
HAZNOW = VAL*pdf1/(1-phi(num1/lambda))                    			            
             			   



;;;;;;;;;; Right censoring only ;;;;;;;;;;;;;;;;;

CS=0
IF (CENSORING.EQ.1.AND.DV.EQ.-1.AND.TIME.GT.0) CS=1                                  ; CS=1 for right censored events
IF (CENSORING.EQ.1.AND.TIME.GT.0) Y = (1 - CS)*HAZNOW*SURV + CS*SURV                 ; for right censored events the likelihood of event at this time or greater is the survival function
					                                             ; for non-censored events the likelihood of event at this time is (hazard function)*(survival function)
;;;;;;;;;; End ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;; Right + interval censoring ;;;;;;;;;;;;;;;;;;;;;;

IF (CENSORING.EQ.2.AND.DV.EQ.-1) Y=SURV

IF (CENSORING.EQ.2.AND.DV.EQ.2) SURVZ=SURV
	 
IF (CENSORING.EQ.2.AND.DV.NE.2) SURVZ=SURVZ

IF (CENSORING.EQ.2.AND.DV.EQ.3) Y=SURVZ-SURV

;;;;;;;;;; End ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
				

$THETA
(0,1 )
(0,1 )
(-10000,1,10000)


$OMEGA
0 FIX 
	

$EST MAXEVAL=9999 METHOD=COND LAPLACE NUMERICAL LIKE SLOW NOABORT NOTHETABOUNDTEST NSIG=3 SIGL=9 PRINT=5 MSFO=runCOMPEV1_101.txt

$COV PRINT=E 

$TABLE ID AGE NEUT PRE_TRE MAX_LEG AUC CMAX DBLOCK TIME DV DVID EVID
AGE_CAT NEUT_CAT PRE_TR_C MAX_CAT AUC_CAT CMAX_CAT
SURV CHAZ HAZNOW FORMAT=s1PE18.9 NOPRINT ONEHEADER FILE=runCOMPEV1_101.tab


  
NM-TRAN MESSAGES 
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
             
 (WARNING  3) THERE MAY BE AN ERROR IN THE ABBREVIATED CODE. THE FOLLOWING
 ONE OR MORE RANDOM VARIABLES ARE DEFINED WITH "IF" STATEMENTS THAT DO NOT
 PROVIDE DEFINITIONS FOR BOTH THE "THEN" AND "ELSE" CASES. IF ALL
 CONDITIONS FAIL, THE VALUES OF THESE VARIABLES WILL BE ZERO.
  
   Y

  
License Registered to: BAST Inc. Ltd
Expiration Date:    14 JAN 2018
Current Date:       20 JUN 2017
Days until program expires : 209
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.3.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 Log_normal model. runCOMPEV1_101.
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      436
 NO. OF DATA ITEMS IN DATA SET:  19
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:  10
 MDV DATA ITEM IS DATA ITEM NO.: 19
0INDICES PASSED TO SUBROUTINE PRED:
  12   9   0   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID AGE NEUT PRE_TRE MAX_LEG AUC CMAX DBLOCK TIME DV DVID EVID AGE_CAT NEUT_CAT PRE_TR_C MAX_CAT AUC_CAT CMAX_CAT MDV
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 CHAZ SURV HAZNOW
0FORMAT FOR DATA:
 (2E3.0,E6.0,2E7.0,E16.0,E7.0,E9.0,E8.0,E2.0,2E4.0/E7.0,5E8.0,1F2.0)

 TOT. NO. OF OBS RECS:      236
 TOT. NO. OF INDIVIDUALS:    200
0LENGTH OF THETA:   3
0DEFAULT THETA BOUNDARY TEST OMITTED:   YES
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
 -0.1000E+05     0.1000E+01     0.1000E+05
0INITIAL ESTIMATE OF OMEGA:
 0.0000E+00
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:             YES
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 SLOW GRADIENT METHOD USED:     YES
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
 HEADER:                YES
 FILE TO BE FORWARDED:   NO
 FORMAT:                s1PE18.9
 LFORMAT:
 RFORMAT:
0USER-CHOSEN ITEMS:
 ID AGE NEUT PRE_TRE MAX_LEG AUC CMAX DBLOCK TIME DV DVID EVID AGE_CAT NEUT_CAT PRE_TR_C MAX_CAT AUC_CAT CMAX_CAT SURV
 CHAZ HAZNOW
1DOUBLE PRECISION PREDPP VERSION 7.3.0

 GENERAL NONLINEAR KINETICS MODEL USING LSODA (ADVAN13)
0MODEL SUBROUTINE USER-SUPPLIED - ID NO. 9999
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         CPT1         ON         YES        YES        YES        YES
    2         OUTPUT       OFF        YES        NO         NO         NO
0NRD VALUE(S) FROM SUBROUTINE TOL:   9
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:     12
   TIME DATA ITEM IS DATA ITEM NO.:          9

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
0ERROR SUBROUTINE INDICATES THAT DERIVATIVES OF COMPARTMENT AMOUNTS ARE USED.
0DES SUBROUTINE USES FULL STORAGE MODE.
1


 #TBLN:      1
 #METH: Laplacian Conditional Estimation

 ESTIMATION STEP OMITTED:                 NO  
 ANALYSIS TYPE:                           POPULATION
 SLOW GRADIENT METHOD USED:               YES 
 CONDITIONAL ESTIMATES USED:              YES 
 CENTERED ETA:                            NO  
 EPS-ETA INTERACTION:                     NO  
 LAPLACIAN OBJ. FUNC.:                    YES 
 NUMERICAL 2ND DERIVATIVES:               YES 
 PRED F SET TO A LIKELIHOOD:              YES 
 NO. OF FUNCT. EVALS. ALLOWED:            9999
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES 
 ESTIMATE OUTPUT TO MSF:                  YES 
 ABORT WITH PRED EXIT CODE 1:             NO  
 IND. OBJ. FUNC. VALUES SORTED:           NO  
 NUMERICAL DERIVATIVE 
       FILE REQUEST (NUMDER):             NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP): 0           
 ETA HESSIAN EVALUATION METHOD (ETADER):  0           
 INITIAL ETA FOR MAP ESTIMATION (MCETA):  0           
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):    9           
 GRADIENT SIGDIGITS OF 
       FIXED EFFECTS PARAMETERS (SIGL):   9           
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   365690.783133644        NO. OF FUNC. EVALS.:   4
 CUMULATIVE NO. OF FUNC. EVALS.:        4
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:  -9.7049E+03 -1.9526E+03  9.9707E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:   360311.149721735        NO. OF FUNC. EVALS.:   6
 CUMULATIVE NO. OF FUNC. EVALS.:       38
 NPARAMETR:  6.7140E-01  2.9572E+02 -5.4830E+00
 PARAMETER: -2.9839E-01  5.7894E+00  9.8703E-02
 GRADIENT:  -1.0183E+02 -6.2753E+01 -2.6907E+03

0ITERATION NO.:   10    OBJECTIVE VALUE:   360254.014263671        NO. OF FUNC. EVALS.:   5
 CUMULATIVE NO. OF FUNC. EVALS.:       64
 NPARAMETR:  1.4412E+00  9.3322E+02  4.3640E+01
 PARAMETER:  4.6548E-01  6.9386E+00  1.0853E-01
 GRADIENT:   5.8763E+00 -9.3822E-01 -3.2834E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:   360253.772970735        NO. OF FUNC. EVALS.:   5
 CUMULATIVE NO. OF FUNC. EVALS.:       92
 NPARAMETR:  1.3995E+00  9.2371E+02  4.6521E+01
 PARAMETER:  4.3609E-01  6.9284E+00  1.0910E-01
 GRADIENT:   7.5190E-02  2.4510E-01 -1.8367E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:   360253.678518176        NO. OF FUNC. EVALS.:   8
 CUMULATIVE NO. OF FUNC. EVALS.:      126
 NPARAMETR:  1.4171E+00  9.8033E+02  5.0648E+01
 PARAMETER:  4.4860E-01  6.9879E+00  1.0993E-01
 GRADIENT:   2.0859E-01 -1.1940E-01 -5.7484E+00

0ITERATION NO.:   22    OBJECTIVE VALUE:   360253.678181308        NO. OF FUNC. EVALS.:   0
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  1.4178E+00  9.8448E+02  5.0878E+01
 PARAMETER:  4.4912E-01  6.9921E+00  1.0998E-01
 GRADIENT:   1.3193E-03 -2.1149E-03  1.0066E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      134
 NO. OF SIG. DIGITS IN FINAL EST.:  4.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         0.0000E+00
 SE:             0.0000E+00
 N:                     200

 P VAL.:         1.0000E+00

 ETAshrink(%):   1.0000E+02
 EBVshrink(%):   0.0000E+00
 EPSshrink(%):   1.0000E-10

 #TERE:
 Elapsed estimation time in seconds:     4.30
 Elapsed covariance time in seconds:     0.64
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************   360253.678       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3     
 
         1.42E+00  9.84E+02  5.09E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


            ETA1   
 
 ETA1
+        0.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


            ETA1   
 
 ETA1
+        0.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3     
 
         1.86E-01  3.10E+02  1.38E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


            ETA1   
 
 ETA1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


            ETA1   
 
 ETA1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11  
 
 TH 1
+        3.46E-02
 
 TH 2
+        4.70E+01  9.62E+04
 
 TH 3
+        5.16E-01  2.27E+03  1.91E+02
 
 OM11
+       ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11  
 
 TH 1
+        1.86E-01
 
 TH 2
+        8.15E-01  3.10E+02
 
 TH 3
+        2.01E-01  5.30E-01  1.38E+01
 
 OM11
+       ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11  
 
 TH 1
+        1.10E+02
 
 TH 2
+       -6.52E-02  5.29E-05
 
 TH 3
+        4.77E-01 -4.53E-04  9.34E-03
 
 OM11
+       ......... ......... ......... .........
 
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************                         LAPLACIAN CONDITIONAL ESTIMATION                       ********************
 ********************                      EIGENVALUES OF COR MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

             1         2         3
 
         1.11E-01  8.18E-01  2.07E+00
 
 #CPUT: Total CPU Time in Seconds,        4.711
Stop Time: 
20/06/2017 
15:01
