$SIZES LIM4=20000 LIM15=20000 
$PROB  Adult Pediatrics Daily and Average Counts Model; Executable_P241
;Negative binomial with Markovian aspects using PDV
;Ahn-type dependency on previous seizures
;Emax CE 
;Eta on Emax
;With Eta on OVDP
;Box-Cox on baseline rates
;Mixture on drug effect
;Change to pediatric subjects: only baseline seizure frequency
;Without subjects with baseline seizures below 2/4 weeks

$DATA  Simulated_P241.csv IGNORE=@
$INPUT STD ID RND VISIT DAYS Q2 CAV DV PDV NDAYS LBAS PED

$CONTR DATA=(LBAS,PED)

$PRED
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
     LTRTE = LPLAC
   ELSE
     LTRTE = LPLAC+LEFF
   ENDIF

   EST     = MIXEST

   LE      = LS0+Q2*LTRTE
   LAMB    = EXP(LE)*NDAYS
   OVDP    = EXP(TVOVDP+PED*ETA(5))
   IF (OVDP.LT.0.0001) OVDP=0.0001
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
 IF (DV.GT.0) THEN
     LDVFAC=(DV+0.5)*LOG(DV)-DV+0.5*LOG(6.283185)+LOG(1+(1/(12*DV)))
 ELSE
     LDVFAC=0
 ENDIF

;Negative binomial:
   Y = -2*(LGAM2-LGAM1-LDVFAC-TRM1+TRM2)

$MIX
   NSPOP=2
   P1PHI=LOG(THETA(9)/(1-THETA(9)))

   COV1=PED*THETA(10)
   P1=EXP(P1PHI+COV1)/(1+EXP(P1PHI+COV1))

   P(1) = P1
   
P(2) = 1-P1

$THETA -1               ;log E0
$THETA (0,3)            ;ES50
$THETA 1.5              ;logP Smax
$THETA -0.2             ;logP Placebo effect
$THETA -3               ;logP Emax 
$THETA  4               ;log EC50 (mg/L)
$THETA -2               ;log Overdispersion Alpha
$THETA  0.5             ;Box-Cox on baseline seizure rate
$THETA (0.01,0.4,0.99)  ;mixture fraction
$THETA 0 FIXED          ;Peds on P1
$THETA 0.4              ;logP Peds on E0
$THETA 0 FIXED          ;logP Peds on Placebo effect
$THETA 0 FIXED          ;logP Peds on Emax
$THETA 0 FIXED          ;logP Peds on EC50

$OMEGA 0.8 1.5 0.2 0.7 10

$EST PRINT=10 MAXEVALS=9999 METHOD=COND NOABORT LAPLACE -2LL
$COV PRINT=E UNCONDITIONAL 
$TABLE STD ID RND VISIT DAYS Q2 CAV DV PDV NDAYS LBAS PED 
       OVDP LS00 ES50 LSMAX LPLAC LEMAX LAMB 
       ETA1 ETA2 ETA3 ETA4 ETA5 EST
       ONEHEADER NOPRINT NOAPPEND FILE=runP241.csv
