$SIZES MAXIDS=20000 LIM1=1200000 LIM2=1200000 LIM3=20000 LIM4=20000 LIM5=20000 LIM6=1200000 LIM7=210 LIM8=6 LIM13=20000 LIM15=20000 LIM16=1200000

$PROB  Adult Pediatrics Daily and Average Counts Model; Simulate_P241
;Negative binomial with Markovian aspects using PDV
;Ahn-type dependency on previous seizures
;Emax CE 
;Eta on Emax
;With Eta on OVDP
;Box-Cox on baseline rates
;Mixture on drug effect
;Change to pediatric subjects: only baseline seizure frequency
;Simulate using Original data structure
;Without subjects with baseline seizures below 2/4 weeks

$DATA  RawData.csv IGNORE=@
$INPUT 
 STD ID VISIT DAYS DV NDAYS DOSE RND Q2 PDV1 PDV=DROP AGE WT HT SEX 
 PED BASE LBAS CL CAV FLGDV7 FLGDV8 FLGDV9 FLGDV10
$DATA IGNORE=(FLGDV7.NE.0)
$DATA IGNORE=(FLGDV8.NE.0)
$DATA IGNORE=(FLGDV9.NE.0)
$DATA IGNORE=(FLGDV10.NE.0)

$CONTR DATA=(PED)

$PRED

;--------Calculate PDV ------------------------------------------------------
;--------PDV is only used for kids (PED=1) with daily seizure counts---------

 IF(NEWIND.NE.2) PREV= 0      ; previous obs at start is set to 0
                 PDV = PREV
 IF (PED.EQ.0)   PDV=-99

;-------- Model -------------------------------------------------------------

   TVLS00  = THETA(1)
   TVES50  = THETA(2)
   TVLSMAX = THETA(3)
   TVLPLAC = THETA(4)
   TVLEMAX = THETA(5)
   TVLEC50 = THETA(6)
   TVOVDP  = THETA(7)
   SHP1    = THETA(8)
   LS00P   = THETA(11)
   LPLACP  = THETA(12)
   LEMAXP  = THETA(13)
   LEC50P  = THETA(14)

;Box-Cox baseline eta
   TETA1   = (EXP(ETA(1))**SHP1-1)/SHP1

   LS00    = TVLS00+TETA1+PED*LS00P
   ES50    = TVES50
   LSMAX   = TVLSMAX+ETA(2)
   LS0     = LS00+PED*LSMAX*PDV/(ES50+PDV)
   LPLAC   = TVLPLAC+ETA(3)+PED*LPLACP
   LEMAX   = TVLEMAX*EXP(ETA(4))+PED*LEMAXP
   LEC50   = TVLEC50+PED*LEC50P
   LEFF    = LEMAX*CAV/(EXP(LEC50)+CAV)

   IF (MIXNUM.EQ.2) THEN
     LTRTE   = LPLAC
   ELSE
     LTRTE   = LPLAC+LEFF
   ENDIF

   LE      = LS0+Q2*LTRTE
   LAMB    = EXP(LE)*NDAYS
   OVDP    = EXP(TVOVDP+PED*ETA(5))
   IF (OVDP.LE.0.0001) OVDP=0.0001
   AGM1    = (OVDP+1)/OVDP
   AGM2    =  DV+1+(1/OVDP)
   PREC1   = (1+(1/(12*AGM1)))
   PREC2   = (1+(1/(12*AGM2)))
   LFAC1   = (AGM1+0.5)*LOG(AGM1)-AGM1+0.5*LOG(6.283185)+PREC1
   LFAC2   = (AGM2+0.5)*LOG(AGM2)-AGM2+0.5*LOG(6.283185)+PREC2
   LGAM1   =  LOG(OVDP)+LOG(OVDP/(1+OVDP))+LFAC1
   LGAM2   = -LOG(DV+(1/OVDP))-LOG(DV+1+(1/OVDP))+LFAC2
   TRM1    = (1/OVDP)*LOG(1+(OVDP*LAMB))
   TRM2    =  DV*LOG(LAMB/(LAMB+(1/OVDP)))

; MODIFIED STIRLING'S APPROXIMATION
 EST=MIXNUM
 IF (DV.GT.0) THEN
     LDVFAC=(DV+0.5)*LOG(DV)-DV+0.5*LOG(6.283185)+LOG(1+(1/(12*DV)))
 ELSE
     LDVFAC=0
 ENDIF

;Negative binomial:
   Y = -2*(LGAM2-LGAM1-LDVFAC-TRM1+TRM2)

;Simulated data set count
 SIMID=IREP

;-------- Simulate data -----------------------------------------------------

 IF (ICALL.EQ.4) THEN
      N=0
      T=0
      ITIMES=0
      CALL RANDOM (2,R)
      DO WHILE (R.GT.T.AND.ITIMES<2400)
        ITIMES=ITIMES+1
        AGM1    = (OVDP+1)/OVDP
        AGM2    =  N+1+(1/OVDP)
        PREC1   = (1+(1/(12*AGM1)))
        PREC2   = (1+(1/(12*AGM2)))
        LFAC1   = (AGM1+0.5)*LOG(AGM1)-AGM1+0.5*LOG(6.283185)+PREC1
        LFAC2   = (AGM2+0.5)*LOG(AGM2)-AGM2+0.5*LOG(6.283185)+PREC2
        LGAM1   =  LOG(OVDP)+LOG(OVDP/(1+OVDP))+LFAC1
        LGAM2   = -LOG(N+(1/OVDP))-LOG(N+1+(1/OVDP))+LFAC2
        TRM1    = (1/OVDP)*LOG(1+(OVDP*LAMB))
        TRM2    =  N*LOG(LAMB/(LAMB+(1/OVDP)))
        IF (N.GT.0) THEN
          LDVFAC=(N+0.5)*LOG(N)-N+0.5*LOG(6.283185)+LOG(1+(1/(12*N)))
        ELSE
          LDVFAC=0
        ENDIF
        YY = EXP(LGAM2-LGAM1-LDVFAC-TRM1+TRM2)
        T=YY+T
        IF (R.GT.T) N=N+1
      END DO
      DV=N
 ENDIF

;-------- Create PDV  for simulated data ------------------------------------
 IF(NEWIND.NE.2) PREV = 0
                 PDV  = PREV
                 PREV = DV     ; saves the current DV value for the next record 
                               ; (which will then be a PDV)
 IF (PED.EQ.0)   PDV=-99

$MIX
   NSPOP=2
   P1PHI=LOG(THETA(9)/(1-THETA(9)))

   COV1=PED*THETA(10)
   P1=EXP(P1PHI+COV1)/(1+EXP(P1PHI+COV1))

   P(1) = P1
   P(2) = 1-P1

$THETA  -1.09E+00  2.75E+00  1.28E+00 -1.60E-01 -3.13E+00  3.45E+00 
        -2.24E+00  4.42E-01  3.35E-01  0 FIXED   4.20E-01  
         0 FIXED 0 FIXED 0 FIXED 

$OMEGA   7.55E-01  1.43E+00  1.66E-01  6.40E-01  8.47E+00

$SIML (987987) (678910 UNI) ONLYSIMULATION NOPREDICTION 

$TABLE STD ID RND VISIT DAYS Q2 CAV DV PDV NDAYS LBAS PED
       NOHEADER NOPRINT NOAPPEND FILE=Simulated_P241.csv
