$SIZES DIMNEW=2000 MAXIDS=5000

$PROBLEM Morphine PD in ventillated newborns using IRT
;        Välitalo et al. CPT:PSP DOI: 10.1002/psp4.12156

$INPUT ID   
       TIME    ; minutes
       DV      ; item score
       ITEM    ; 1. COMFORT alertness
               ; 2. COMFORT calmness/agitation
               ; 3. COMFORT respiratory response
               ; 5. COMFORT body movement
               ; 7. COMFORT facial tension
               ; 12. VAS score (cm, range 0 - 10)
               ; 25. PIPP brow bulge
               ; 26. PIPP eye squeeze
               ; 27. PIPP nasolabial furrow
               ; 28. NIPS
       CP      ; morphine plasma concentration
       MDV
       MOMENT  ; 1. before endotrachial succtioning
               ; 2. during endotrachial succtioning
               ; 3. after endotrachial succtioning
       OBSTYPE ; 1. observation based on video
               ; 2. bedside observation by nurse

$ABBREVIATED DECLARE INTEGER IND 

$DATA Simulated_MorphineVentillatedNewborns.csv IGNORE=@

$PRED  
EFF=THETA(56)*CP            ; morphine effect
COVEFF=THETA(57)*TIME/24/60 ; time effect

; identify time of pain measurement relative to endotracheal succtioning (pre / during / after)
F_FLAG=1  
MU_1=THETA(53)
presuctpain=MU_1+ETA(1)
MU_2=THETA(54)
suctpain=MU_2+ETA(2)
MU_3=THETA(55)
aftsuctpain=MU_3+ETA(3)

PAIN=presuctpain
IF (moment.EQ.2) PAIN=suctpain
IF (moment.EQ.3) PAIN=aftsuctpain
PAIN=PAIN-EFF+COVEFF

TVPAIN=THETA(53)
IF (moment.EQ.2) TVPAIN=THETA(54)
IF (moment.EQ.3) TVPAIN=THETA(55)
TVPAIN=TVPAIN-EFF+COVEFF

IND=1
IF(ITEM.EQ.2) IND=6
IF(ITEM.EQ.3) IND=11
IF(ITEM.EQ.5) IND=16
IF(ITEM.EQ.7) IND=21
IF(ITEM.EQ.25) IND=26
IF(ITEM.EQ.26) IND=31
IF(ITEM.EQ.27) IND=36
IF(ITEM.EQ.28) IND=41

discr=THETA(IND)   
diff1=THETA(IND+1) 
diff2=THETA(IND+2) 
diff3=THETA(IND+3) 
diff4=THETA(IND+4) 
diff5=THETA(IND+5)
diff6=THETA(IND+6) 
diff7=THETA(IND+7)

diffGT1=diff1
diffGT2=diff1+diff2
diffGT3=diff1+diff2+diff3
diffGT4=diff1+diff2+diff3+diff4
diffGT5=diff1+diff2+diff3+diff4+diff5
diffGT6=diff1+diff2+diff3+diff4+diff5+diff6
diffGT7=diff1+diff2+diff3+diff4+diff5+diff6+diff7

pGT1=EXP(discr*(PAIN-diffGT1))/(1+EXP(discr*(PAIN-diffGT1)))
pGT2=EXP(discr*(PAIN-diffGT2))/(1+EXP(discr*(PAIN-diffGT2)))
pGT3=EXP(discr*(PAIN-diffGT3))/(1+EXP(discr*(PAIN-diffGT3)))
pGT4=EXP(discr*(PAIN-diffGT4))/(1+EXP(discr*(PAIN-diffGT4)))
pGT5=EXP(discr*(PAIN-diffGT5))/(1+EXP(discr*(PAIN-diffGT5)))
pGT6=EXP(discr*(PAIN-diffGT6))/(1+EXP(discr*(PAIN-diffGT6)))
pGT7=EXP(discr*(PAIN-diffGT7))/(1+EXP(discr*(PAIN-diffGT7)))

peq1=1-pGT1
peq2=pGT1-pGT2
peq3=pGT2-pGT3
peq4=pGT3-pGT4
peq5=pGT4

IF (ITEM.EQ.28) THEN ; NIPS
peq5=pGT4-pGT5
peq6=pGT5-pGT6
peq7=pGT6-pGT7
peq8=pGT7
ENDIF

IF(peq1.LT.0) peq1=0
IF(peq2.LT.0) peq2=0
IF(peq3.LT.0) peq3=0
IF(peq4.LT.0) peq4=0
IF(peq5.LT.0) peq5=0
IF(peq6.LT.0) peq6=0
IF(peq7.LT.0) peq7=0
IF(peq8.LT.0) peq8=0

COMFADD=0
IF (ITEM.LT.10) COMFADD=1 
IF(DV.EQ.0+COMFADD) P=peq1
IF(DV.EQ.1+COMFADD) P=peq2
IF(DV.EQ.2+COMFADD) P=peq3
IF(DV.EQ.3+COMFADD) P=peq4
IF(DV.EQ.4+COMFADD) P=peq5
IF(DV.EQ.5+COMFADD) P=peq6
IF(DV.EQ.6+COMFADD) P=peq7
IF(DV.EQ.7+COMFADD) P=peq8

Y=P

IF (ITEM.LT.10) THEN
IEXPECT=1*peq1+2*peq2+3*peq3+4*peq4+5*peq5
IEVAR=(1-IEXPECT)**2*peq1+ (2-IEXPECT)**2*peq2+ (3-IEXPECT)**2*peq3+ (4-IEXPECT)**2*peq4+(5-IEXPECT)**2*peq5
ELSE
IEXPECT=0*peq1+1*peq2+2*peq3+3*peq4+4*peq5+5*peq6+6*peq7+7*peq8
IEVAR1=(0-IEXPECT)**2*peq1+ (1-IEXPECT)**2*peq2+ (2-IEXPECT)**2*peq3+ (3-IEXPECT)**2*peq4
IEVAR=IEVAR1+ (4-IEXPECT)**2*peq5+ (5-IEXPECT)**2*peq6+ (6-IEXPECT)**2*peq7+ (7-IEXPECT)**2*peq8
ENDIF

;tvpredicts
tpGT1=EXP(discr*(TVPAIN-diffGT1))/(1+EXP(discr*(TVPAIN-diffGT1)))
tpGT2=EXP(discr*(TVPAIN-diffGT2))/(1+EXP(discr*(TVPAIN-diffGT2)))
tpGT3=EXP(discr*(TVPAIN-diffGT3))/(1+EXP(discr*(TVPAIN-diffGT3)))
tpGT4=EXP(discr*(TVPAIN-diffGT4))/(1+EXP(discr*(TVPAIN-diffGT4)))
tpGT5=EXP(discr*(TVPAIN-diffGT5))/(1+EXP(discr*(TVPAIN-diffGT5)))
tpGT6=EXP(discr*(TVPAIN-diffGT6))/(1+EXP(discr*(TVPAIN-diffGT6)))
tpGT7=EXP(discr*(TVPAIN-diffGT7))/(1+EXP(discr*(TVPAIN-diffGT7)))
tpeq1=1-tpGT1
tpeq2=tpGT1-tpGT2
tpeq3=tpGT2-tpGT3
tpeq4=tpGT3-tpGT4
tpeq5=tpGT4
tpeq6=0
tpeq7=0
tpeq8=0

IF (ITEM.EQ.28) THEN ; NIPS
tpeq5=tpGT4-tpGT5
tpeq6=tpGT5-tpGT6
tpeq7=tpGT6-tpGT7
tpeq8=tpGT7
ENDIF

IF (ITEM.LT.10) THEN
EXPECT=1*tpeq1+2*tpeq2+3*tpeq3+4*tpeq4+5*tpeq5
EVAR=(1-EXPECT)**2*tpeq1+ (2-EXPECT)**2*tpeq2+ (3-EXPECT)**2*tpeq3+ (4-EXPECT)**2*tpeq4+(5-EXPECT)**2*tpeq5
ELSE
EXPECT=0*tpeq1+1*tpeq2+2*tpeq3+3*tpeq4+4*tpeq5+5*tpeq6+6*tpeq7+7*tpeq8
EVAR1=(0-EXPECT)**2*tpeq1+ (1-EXPECT)**2*tpeq2+ (2-EXPECT)**2*tpeq3+ (3-EXPECT)**2*tpeq4
EVAR=EVAR1+ (4-EXPECT)**2*tpeq5+ (5-EXPECT)**2*tpeq6+ (6-EXPECT)**2*tpeq7+ (7-EXPECT)**2*tpeq8
ENDIF

VASDIFF=THETA(49)
IF (OBSTYPE.EQ.2)  VASDIFF=THETA(50) ;bedside
VASDISCR=THETA(51)
IF (OBSTYPE.EQ.2)  VASDISCR=THETA(52) ;bedside
VASPRED=10*EXP(VASDISCR*(PAIN-VASDIFF))/(1+EXP(VASDISCR*(PAIN-VASDIFF)))
VASTYPPRED=10*EXP(VASDISCR*(TVPAIN-VASDIFF))/(1+EXP(VASDISCR*(TVPAIN-VASDIFF)))

IF (ITEM.EQ.12) F_FLAG=0
IF (ITEM.EQ.12) Y=VASPRED+EPS(1)*(2-OBSTYPE)+EPS(2)*(OBSTYPE-1)

$ESTIMATION METHOD=1 PRINT=1 MAXEVAL=9999 LAPLACE
$COVARIANCE UNCONDITIONAL PRINT=E
$THETA  
 (0,3.6) FIX   ; discrimination COMFORT alertness
 -0.761 FIX    ; 11.   lowest difficulty item 1 
 (0,1.39) FIX  ; 12.   second lowest difficulty item 1 
 (0,0.887) FIX ; 13.   second highest difficulty item 1 
 (0,1.01) FIX  ; 14.   highest difficulty item 1
 (0,4.43) FIX  ; 2.discrimination COMFORT calmness/agitation,
 -0.162 FIX    ; 15.   lowest difficulty item 2 
 (0,0.992) FIX ; 16.   second lowest difficulty item 2 
 (0,0.808) FIX ; 17.   second highest difficulty item 2
 (0,0.533) FIX ; 18.   highest difficulty item 2
 (0,2.91) FIX  ; 3.discrimination COMFORT respiratory response,
 -2.67 FIX     ; 19.   lowest difficulty item 3 
 (0,3.2) FIX   ; 20.   second lowest difficulty item 3
 (0,1.72) FIX  ; 21.   second highest difficulty item 3
 1.09 FIX      ; 22.   highest difficulty item 3
 (0,2.89) FIX  ; 5.discrimination COMFORT body movement,
 -1.74 FIX     ; 27.   lowest difficulty item 5 
 (0,0.777) FIX ; 28.   second lowest difficulty item 5
 (0,1.94) FIX  ; 29.   second highest difficulty item 5
 (0,1.94) FIX  ; 30.   highest difficulty item 5
 (0,4.04) FIX  ; 7.discrimination COMFORT facial tension,
 -1.49 FIX     ; 35.   lowest difficulty item 7 
 (0,1.76) FIX  ; 36.   second lowest difficulty item 7
 (0,1.43) FIX  ; 37.   second highest difficulty item 7
 (0,0.334) FIX ; 38.   highest difficulty item 7

 (0,2.21) FIX  ; discrimination pipp brow bulge
 0.467 FIX     ; lowest difficulty item 25 
 (0,1.08) FIX  ; second lowest difficulty item 25
 (0,0.937) FIX ; highest difficulty item 25
 100000 FIX    ; difficulties pipp5_
 (0,2.14) FIX  ; discrimination pipp eye squeeze
 0.343 FIX     ; lowest difficulty item 26 
 (0,1.17) FIX  ; second lowest difficulty item 26
 (0,0.98) FIX  ; highest difficulty item 26
 100000 FIX    ; difficulties pipp6_
 (0,2.66) FIX  ; discrimination pipp nasolabial furrow
 0.442 FIX     ; lowest difficulty item 27 
 (0,0.829) FIX ; second lowest difficulty item 27
 (0,1.27) FIX  ; highest difficulty item 27
 100000 FIX    ; difficulties pipp7_

 (0,4.37) FIX  ; discrimination nips
 -0.082 FIX    ; lowest difficulty item 28 
 (0,0.332) FIX ; difficulty item 28
 (0,0.239) FIX ; difficulty item 28
 (0,0.216) FIX ; difficulty item 28
 (0,0.211) FIX ; difficulty item 28
 (0,0.249) FIX ; difficulty item 28
 (0,0.492) FIX ; higest difficulty item 28

 1.907720 FIX  ; DIFFICULTY VAS, ETA VALUE TO GIVE 50 % RESPONSE
 2.576450 FIX  ; VAS difficulty bedside
 (0,1.160130) FIX ; VAS DISCRIMINATION
 (0,0.663043) FIX ; VAS discrimination bedside

 -0.27         ; presuctioning
 1.15          ; suctioning
 -0.393        ; aftersuctioning
 0.0091        ; morphine effect slope
 (0,0.0476)    ; increase in pain by time

$OMEGA  BLOCK(3)
 0.311  
 0.172 0.34
 0.216 0.153 0.317

$SIGMA  0.755
$SIGMA  2.55

$TABLE  ID TIME DV ITEM CP MDV MOMENT OBSTYPE NOPRINT ONEHEADER FILE=sdtab1

