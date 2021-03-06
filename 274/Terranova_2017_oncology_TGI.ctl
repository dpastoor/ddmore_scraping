; Script generated by the pharmML2Nmtran Converter v.0.4.0
; Source 	: PharmML 0.8.1
; Target 	: NMTRAN 7.3.0
; Model 	: Terranova_2017_oncology_TGI
; Model Name: Generated from MDL. MOG ID: DEB_TGI_mog
; Dated 	: Thu Dec 21 15:21:36 CET 2017

$PROBLEM Generated from MDL. MOG ID: DEB_TGI_mog

$INPUT  ID TIME DV DVID AMT EVID CMT
$DATA "Simulated_DEB_TGI_data.csv" IGNORE=@
$SUBS ADVAN13 TOL=9

$MODEL 
COMP (COMP1) 	;Z
COMP (COMP2) 	;Q1
COMP (COMP3) 	;Q2
COMP (COMP4) 	;EN
COMP (COMP5) 	;VU1
COMP (COMP6) 	;VU2
COMP (COMP7) 	;VU3
COMP (COMP8) 	;VU4



$PK 
MU_POP = THETA(1)
MU_U_POP = THETA(2)
GU_POP = THETA(3)
DELTA_VMAX_POP = THETA(4)
W_INITIAL_POP = THETA(5)
VU1_INITIAL_POP = THETA(6)
IC50_POP = THETA(7)
K1_POP = THETA(8)
K2_POP = THETA(9)
B_W = THETA(10)
B_WU = THETA(11)
K10_POP = THETA(12)
K12_POP = THETA(13)
K21_POP = THETA(14)
V1_POP = THETA(15)
EN_INITIAL_POP = THETA(16)
XI_POP = THETA(17)
NI_POP = THETA(18)
GR_POP = THETA(19)
V1INF_POP = THETA(20)
RHO_B_POP = THETA(21)



K10 = K10_POP

K12 = K12_POP

K21 = K21_POP

V1 = V1_POP

EN_INITIAL = EN_INITIAL_POP

RHO_B = RHO_B_POP

XI = XI_POP

NI = NI_POP

GR = GR_POP

V1INF = V1INF_POP

MU = MU_POP

MU_U = MU_U_POP

GU = GU_POP

DELTA_VMAX = DELTA_VMAX_POP

W_INITIAL = W_INITIAL_POP

VU1_INITIAL = VU1_INITIAL_POP

IC50 = IC50_POP

K1 = K1_POP

K2 = K2_POP

DENSITY_V = 1

DENSITY_VU = 1

OMEG = 0.75

M = (NI/((V1INF**(1/3))*GR))

Z_INITIAL = (W_INITIAL/(1+(EN_INITIAL*XI)))

A_0(1) = Z_INITIAL
A_0(2) = 0
A_0(3) = 0
A_0(4) = EN_INITIAL
A_0(5) = VU1_INITIAL
A_0(6) = 0
A_0(7) = 0
A_0(8) = 0

$DES 
Z_DES = A(1)
Q1_DES = A(2)
Q2_DES = A(3)
EN_DES = A(4)
VU1_DES = A(5)
VU2_DES = A(6)
VU3_DES = A(7)
VU4_DES = A(8)
C_DES = (Q1_DES/V1)
RHO_DES = (RHO_B*(1-(C_DES/(IC50+C_DES))))
KU_DES = ((MU_U*VU1_DES)/(Z_DES+(MU_U*VU1_DES)))
SWITCH1_DES = ((((((1-KU_DES)*NI)*EN_DES)*(Z_DES**(2/3)))-((GR*M)*Z_DES))/(GR+((1-((MU_U*VU1_DES)/(Z_DES+(MU_U*VU1_DES))))*EN_DES)))
SWITCH2_DES = ((((((1-KU_DES)*NI)*EN_DES)*(Z_DES**(2/3)))-((GR*M)*Z_DES))/((1-((MU_U*VU1_DES)/(Z_DES+(MU_U*VU1_DES))))*(EN_DES+(OMEG*GR))))
WU_DES = (DENSITY_VU*(((VU1_DES+VU2_DES)+VU3_DES)+VU4_DES))
W_DES = ((DENSITY_V*(1+(XI*EN_DES)))*Z_DES)
W_ERR_DES = (B_W*SQRT(W_DES))
WU_ERR_DES = (B_WU*SQRT(WU_DES))
DADT(2) = ((K21*Q2_DES)-((K10+K12)*Q1_DES))
DADT(3) = ((K12*Q1_DES)-(K21*Q2_DES))
DADT(4) = ((NI/(Z_DES**(1/3)))*((RHO_DES*((V1INF/(VU1_DES+Z_DES))**(2/3)))-EN_DES))
DADT(5) = ((((KU_DES/GU)*((((EN_DES*NI)*(Z_DES**(2/3)))+(DELTA_VMAX*EN_DES))+((DELTA_VMAX*OMEG)*GR)))-(MU*VU1_DES))-((K2*C_DES)*VU1_DES))
DADT(6) = (((K2*C_DES)*VU1_DES)-(K1*VU2_DES))
DADT(7) = ((K1*VU2_DES)-(K1*VU3_DES))
DADT(8) = ((K1*VU3_DES)-(K1*VU4_DES))

$ERROR 
EPS_RES_W = EPS(1)
EPS_RES_WU = EPS(2)

Z = A(1)
Q1 = A(2)
Q2 = A(3)
EN = A(4)
VU1 = A(5)
VU2 = A(6)
VU3 = A(7)
VU4 = A(8)
C = (Q1/V1)
RHO = (RHO_B*(1-(C/(IC50+C))))
KU = ((MU_U*VU1)/(Z+(MU_U*VU1)))
SWITCH1 = ((((((1-KU)*NI)*EN)*(Z**(2/3)))-((GR*M)*Z))/(GR+((1-((MU_U*VU1)/(Z+(MU_U*VU1))))*EN)))
SWITCH2 = ((((((1-KU)*NI)*EN)*(Z**(2/3)))-((GR*M)*Z))/((1-((MU_U*VU1)/(Z+(MU_U*VU1))))*(EN+(OMEG*GR))))
WU = (DENSITY_VU*(((VU1+VU2)+VU3)+VU4))
W = ((DENSITY_V*(1+(XI*EN)))*Z)
W_ERR = (B_W*SQRT(W))
WU_ERR = (B_WU*SQRT(WU))

IF(DVID.EQ.1) THEN
	IPRED = W
W = W_ERR
Y = IPRED+W*EPS_RES_W
IRES = DV - IPRED
IWRES = IRES/W

ENDIF

IF(DVID.EQ.2) THEN
	IPRED = WU
W = WU_ERR
Y = IPRED+W*EPS_RES_WU
IRES = DV - IPRED
IWRES = IRES/W

ENDIF

$THETA 
( 0.0 , 0.0223 )	;mu_POP
( 0.0 , 13.3 )	;mu_u_POP
( 0.0 , 11.7 )	;gu_POP
( 0.0 , 0.185 )	;delta_Vmax_POP
( 0.0 , 21.2 )	;W_initial_POP
( 0.0 , 0.0023 )	;Vu1_initial_POP
( 0.0 , 0.461 )	;IC50_POP
( 0.0 , 0.462 )	;k1_POP
( 0.0 , 6.53E-4 )	;k2_POP
( 0.0 , 0.101 )	;b_W
( 0.0 , 0.134 )	;b_Wu
(20.832  FIX )	;K10_POP
(0.144  FIX )	;K12_POP
(2.011  FIX )	;K21_POP
(813.1  FIX )	;V1_POP
(1.3  FIX )	;En_initial_POP
(0.184  FIX )	;xi_POP
(1.2242  FIX )	;ni_POP
(12.2  FIX )	;gr_POP
(22.6  FIX )	;V1inf_POP
(1.0  FIX )	;rho_b_POP

$OMEGA 0 FIX

$SIGMA 
1.0 FIX
1.0 FIX


$EST METHOD=COND NSIG=3 SIGL=9 MAXEVALS=9999 PRINT=10 NOABORT

$TABLE  ID TIME DVID AMT EVID CMT PRED RES WRES DV IPRED IRES IWRES Y NOAPPEND NOPRINT FILE=sdtab

$TABLE  ID K10 K12 K21 V1 EN_INITIAL RHO_B XI NI GR V1INF MU MU_U GU DELTA_VMAX W_INITIAL VU1_INITIAL IC50 K1 K2 DENSITY_V DENSITY_VU OMEG M Z_INITIAL NOAPPEND NOPRINT FILE=patab


