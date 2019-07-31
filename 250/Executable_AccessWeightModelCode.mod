$PROBLEM Mida PK in obese adolescents and adults with access weight quantification
$INPUT ID   ;<=30 = adults, >30 = adolescent
       TIME ;in min
       AMT  ;in microgram
       RATE ;in microgram / min
       DV   ;in microgram / L
       MDV  
       CMT 
       TBW  ;total bodyweight in kg
       WTAL ;weight for age and length in kg
       WTAC ;access weight

$DATA FinalDatasetMidaObeseAccessWeight.csv IGNORE=@  
$SUBROUTINES ADVAN6 TOL=5         

$MODEL
COMP=(PODOSE)
COMP=(CENTRAL)
COMP=(PERIP)
COMP=(TRANSIT1)
COMP=(TRANSIT2)
COMP=(TRANSIT3)
COMP=(TRANSIT4)
COMP=(TRANSIT5)

$PK 							
IF(ID.LE.30) TVCL=THETA(1) 
IF(ID.GT.30) TVCL=(THETA(1)*(WTAL/70)**THETA(7))+(THETA(8)*WTAC) 
CL=TVCL*EXP(ETA(1))		;clearance L/min
F1=THETA(2)*EXP(ETA(2))		;bioavailability
V2=THETA(3)*EXP(ETA(5))        	;central volume L 
Q=THETA(4)*EXP(ETA(3))         	;intercompartmental CL L/min 
IF(ID.LE.30) TVV3=THETA(5)*(TBW/141.8)**THETA(9) 
IF(ID.GT.30) TVV3=THETA(5)
V3= TVV3*EXP(ETA(6))            ;peripheral volume L
KA= THETA(6)*EXP(ETA(4))	;absoprtion rate min-1
KTR= KA                         ;transit rate min-1

S2=V2                            
S3=V3

K14=KA
K45=KTR
K56=KTR
K67=KTR
K78=KTR
K82=KTR
K20=CL/V2
K23=Q/V2
K32=Q/V3

$DES
DADT(1)= -K14*A(1)
DADT(2)= KTR*A(8) -K23*A(2) +K32*A(3) -K20*A(2)
DADT(3)= K23*A(2) -K32*A(3)
DADT(4)= K14*A(1) -KTR*A(4)
DADT(5)= KTR*A(4) -KTR*A(5)
DADT(6)= KTR*A(5) -KTR*A(6)
DADT(7)= KTR*A(6) -KTR*A(7)
DADT(8)= KTR*A(7) -KTR*A(8)

$ERROR
IPRED=F                    
Y=F*(1+ERR(1))     ;proportional error model

IRES=DV-IPRED            
DEL=0                         
IF(IPRED.EQ.0)DEL=1
IWRES=(1-DEL)*IRES/(IPRED+DEL)   

$THETA
(0, 0.447) ;CL ADULTS (L/min)
(0, 0.56) ;F1 
(0, 53.7) ;V2 (L)
(0, 1.15) ;Q (L/min)
(0, 168) ;V3 (L)
(0, 0.114) ;KA 
(0.75) FIX ;POW CL
(0, 0.00698) ;overweight 
(0, 3.22) ;TBW POW V3 ADULTS

$OMEGA              
 0.0497 ;CL
 0.143 ;F1
 0.175 ;Q 
 0.226 ;KA
 0.259 ;V2 
 0.146 ;V3 

$SIGMA      
 0.0892 ;porportional

$EST SIGDIG=3 MAXEVAL=9999 PRINT=5 NOABORT METHOD=1 INTERACTION POSTHOC
$COV 
$TABLE ID TIME IPRED IWRES CWRES CMT NOPRINT ONEHEADER NOAPPEND FILE=sdtab143

