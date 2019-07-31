[PROB]
# Zierhut et al. (2008)
  - Population PK-PD model for Fc-osteoprotegerin in healthy postmenopausal women
  - J Pharmacokinet Pharmacodyn 35(4): 379-99 


[PARAM]  @annotated
IV : 0 : IV dose indicator

[CMT] @annotated
SC   : Subcutaneous dosing compartment (mg)
CENT : Central compartment (mg)
P1   : First peripheral compartment (mg)
P2   : Second peripheral compartment (mg)
NTX  : Urinary N-telopeptide
 
[PARAM] @annotated
TVCL   : 168 : Clearance (ml/h)
TVVC   : 2800   : Central volume (ml)
TVVP1  : 443    : Volume of first peripheral cmt (ml)
TVVP2  : 269    : Volume of second peripheral cmt (ml)
TVQ1   : 15.5   : Distribution clearance (ml/h)
TVQ2   : 3.02   : Distribution clearance (ml/h)
TVKA   : 0.0131 : Absorption rate constant (1/h)
TVVMAX : 13300  : Maximum velocity (ng/h)
TVKM   : 6.74   : Michaelis constant (ng/ml)
TVFSC  : 0.0719 : Bioavailability of SC dose (.)
TVKSYN : 0.864  : Biomarker synthesis rate (.)
TVKDEG : 0.0204 : Biomarker elimination rate constant (1/h)
TVIC50 : 5.38   : Half-maximal inhibitory conc (ng/ml)

[GLOBAL]
#define CP (CENT/(VC/1000000.0))

[MAIN]
double CL   = exp(log(TVCL)  + ECL);
double VC   = exp(log(TVVC)  + EVC);
double VP1  = exp(log(TVVP1) + EVP1);
double VP2  = exp(log(TVVP2) + EVP2);
double Q1   = exp(log(TVQ1)  + EQ1);
double Q2   = TVQ2;
double KA   = exp(log(TVKA)  + EKA);
double VMAX = TVVMAX;
double KM   = TVKM;
double FSC  = exp(log(TVFSC) + EFSC);
double KSYN = exp(log(TVKSYN) + EKSYN);
double KDEG = exp(log(TVKDEG) + EKDEG);
double IC50 = exp(log(TVIC50) + EIC50);

NTX_0 = KSYN/KDEG;

F_SC = FSC/(1.0+FSC);

[OMEGA] @annotated
ECL  : 0.0391 : IIV on clearance
EVC  : 0.0102 : IIV on VC
EVP1 : 0.0144 : IIV on VP1
EVP2 : 0.0333 : IIV on VP2
EQ1  : 0.0379 : IIV on Q1
EKA  : 0.0457 : IIV on KA
EFSC : 0.263  : IIV on FSC

[OMEGA] @block @annotated
EKSYN : 0.281               : IIV on KSYN
EKDEG : 0.0867 0.0325       : IIV on KDEG
EIC50 : 0.0000 0.0000  1.18 : IIV on IC50
  
[SIGMA] @annotated
ADDIV : 0.0193 : Additive error IV dose
ADDSC : 0.7330 : Additive error SC dose
  
[SIGMA] @annotated
PDPROP : 0.0407 : Proportional error NTX
PDADD  :   20.7 : Additive error NTX

[ODE]
double CLNL = VMAX/(CP+KM);
dxdt_SC     = -KA*SC;
dxdt_CENT   =  KA*SC - (CL+Q1+Q2+CLNL)*CENT/VC + Q1*P1/VP1 + Q2*P2/VP2;
dxdt_P1     =  CENT*Q1/VC - P1*Q1/VP1;
dxdt_P2     =  CENT*Q2/VC - P2*Q2/VP2;
dxdt_NTX    =  KSYN*(1.0 - CP/(IC50+CP)) - KDEG*NTX;

[TABLE]
double IPRED = CP;
double PKEPS = IV==1 ? ADDIV : ADDSC;
double PKDV = exp(log(IPRED)+PKEPS);
double PDDV = NTX*(1+PDPROP) + PDADD;

[CAPTURE] @annotated
PKDV : Fc-OPG serum concentration (ng/ml)
PDDV : NTX (nM BCE/mM creatinine per hr)

  
