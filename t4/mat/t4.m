close all
clear all
format long

f = 100e3;
w = 2*pi*f;


%gain stage

VT=25e-3
BFN=178.7
VAFN=69.7
VBEON=0.7


 printf("ModelData_TAB \n"); 
 printf("$V_{A}$ = %f V\n", VAFN);
 printf("$V_{BEON}$ = %f V \n", VBEON);
 printf("$V_{T}$ = %f V\n", VT);
 printf("$beta$ = %f \n", BFN); 
 printf("ModelData_END \n \n");


RE1=100
RC1=1000
RB1=80000
RB2=20000
VCC=12
RS=100

Ci=1e-3;
Cb=1e-3;
Co=1e-6;
RL=8;
Vinput=1;


 printf("UsedValues_TAB \n"); 
 printf("$V_{CC}$ = %e V\n", VCC);
 printf("$R_{inp}$ = %e Ohm\n", RS); 
 printf("$R_{1}$ = %e Ohm \n", RB1);
 printf("$R_{2}$ = %e Ohm \n", RB2);
 printf("$R_{c}$ = %e Ohm \n", RC1);
 printf("$R_{e}$ = %e Ohm \n", RE1);
 printf("$C_{inp}$ = %e F \n", Ci);
 printf("$C_{b}$ = %e F \n", Cb);
 printf("UsedValues_END \n \n");


%op
RB=1/(1/RB1+1/RB2)
VEQ=RB2/(RB1+RB2)*VCC
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1)
IC1=BFN*IB1
IE1=(1+BFN)*IB1
VE1=RE1*IE1
VO1=VCC-RC1*IC1
VCE=VO1-VE1

 printf("GainOP_TAB \n"); 
 printf("$R_{B}$ = %e Ohm\n", RB);
 printf("$V_{eq}$ = %e V\n", VEQ); 
 printf("$I_{B1}$ = %e A \n", IB1);
 printf("$I_{C1}$ = %e A \n", IC1);
 printf("$I_{E1}$ = %e A \n", IE1);
 printf("$V_{E1}$ = %e V \n", VE1);
 printf("$V_{O1}$ = %e V \n", VO1);
 printf("$V_{CE}$ = %e V \n", VCE);
 printf("GainOP_END \n \n");

%incremental
gm1=IC1/VT
rpi1=BFN/gm1
ro1=VAFN/IC1

RSB=RB*RS/(RB+RS)

%ac
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple = RB/(RB+RS) * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=0
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))

 printf("GainAV_TAB \n"); 
 printf("$A_{V1}$ = %e dB\n", AVI_DB);
 printf("$A_{VIsimple}$ = %e dB\n", AVIsimple_DB); 
 printf("GainAV_END \n \n"); 
 
RE1=100
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)))
ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) ) 
ZO1 = 1/(1/ZX+1/RC1)

RE1=0
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZO1 = 1/(1/ro1+1/RC1)

 printf("GainImpedance_TAB \n"); 
 printf("$Z_{In1}$ = %e Ohm\n", ZI1);
 printf("$Z_{Oot1}$ = %e Ohm\n", ZO1); 
 printf("GainImpedance_END \n \n");
 
%ouput stage
BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

 printf("OutputModel_TAB \n"); 
 printf("$beta$ = %f \n", BFP); 
 printf("$V_{AFP}$ = %f V\n", VAFP);
 printf("$V_{BEON}$ = %f V \n", VEBON);
 printf("OutputModel_END \n \n");
 
 printf("OutputOP_TAB \n"); 
 printf("$V_{I2}$ = %e V\n", VI2);
 printf("$I_{E2}$ = %e A \n", IE2);
 printf("$I_{C2}$ = %e A \n", IC2);
 printf("$V_{O2}$ = %e V \n", VO2);
 printf("OutputOP_END \n \n");
 

gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2

AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)

rpi2 = BFP/gm2;
ro2 = VAFP/IC2;
gm2 = IC2/VT;


 printf("OutputAV_TAB \n"); 
 printf("$A_{V2}$ = %e dB\n", AV2);
 printf("OutputAV_END \n \n"); 
 
 printf("OutputImpedance_TAB \n"); 
 printf("$Z_{In2}$ = %e Ohm\n", ZI2);
 printf("$Z_{Out2}$ = %e Ohm\n", ZO2); 
 printf("OutputImpedance_END \n \n");



%total
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)


 printf("Total_TAB \n"); 
 printf("$A_{V}$ = %e dB\n", AV_DB);
 printf("$Z_{I}$ = %e Ohm\n", ZI);
 printf("$Z_{O}$ = %e Ohm\n", ZO); 
 printf("Total_END \n \n");
 
 printf("ZTOTAL_TAB \n"); 
 printf("$Z_{I}$ = %e Ohm\n", ZI);
 printf("$Z_{O}$ = %e Ohm\n", ZO); 
 printf("ZTOTAL_END \n \n");



% GRANDEZAS PARA PÔR NAS TABELAS

%tabela zin/zout/gain
diary on
Zi = ZI
Zo = ZO
Zi1 = ZI1
Zo1 = ZO1
Zi2 = ZI2
Zo2 = ZO2
Gain = abs(AV)
Gain1 = abs(AV1)
Gain2 = abs(AV2)
diary off

%tabela op
diary on
Vbase = VEQ
Vcoll = -RC1*IC1+VCC
Vemit = VE1
Vemit2 = -RE2*IE2+VCC
Vin = 0
Vin2 = 0
Vout = 0
Vvcc = VCC
diary off

% FREQUENCY RESPONSE
f_H=f;
f_L= (1/(3*Ci)+ 1/((ZO+RL)*Co) +1/((ZI+RS)*Cb))/(2*pi);
band=f_H - f_L

 printf("LC_TAB \n"); 
 printf("$Lower CO freq$ = %e Hz\n", f_L);
 printf("LC_END \n \n");

freq=logspace(1,8);

for i=1:50
wfreq=2*pi*freq(i);
ZCI   = 1/(wfreq*j*Ci);
ZCB   = 1/(wfreq*j*Cb);
ZEB = 1/(1/RE1 + 1/ZCB);
zpi2  = 1/gpi2;
zo2   = 1/go2;
ZE2 = 1/(1/zo2 + 1/RE2);
ZCO   = 1/(wfreq*j*Co);


A = [ RS+ZCI+RB , -RB, 0 , 0 , 0 , 0 , 0; ...
 	    0 , -ZEB , -ro1, ZEB + ro1 + RC1, -RC1, 0, 0 ; ...
      0 , rpi1*gm1 , 1  , 0 , 0, 0 , 0  ; ...
      0 , 0 , 0  , -RC1, zpi2+RC1, ZE2 , -ZE2; ...
      0 , 0  , 0  , 0  , -1-zpi2*gm2 , 1  , 0 ; ...
      0, 0 , 0  , 0 , 0 , -ZE2, ZE2+ZCO+RL ; ...
      RB  , RB + rpi1 + ZEB , 0 , -ZEB  , 0 , 0 , 0 ];
        
B = [Vinput; 0; 0 ; 0 ; 0; 0; 0];

C = A\B;

Gain(i) = abs(C(7)*RL/Vinput);
GaindB(i) = 20*log10(abs(C(7)*RL/Vinput));

endfor

% Plot do gain com a frequência
hfa = figure;
semilogx(freq, Gain,";gain(f);");
xlabel ("f [Hz]");
ylabel ("gain");
print (hfa, "gain_octave.odg", "-depsc");

% Plot do gain em dB com a frequência
hfb = figure;
semilogx(freq, GaindB,";gaindB(f);");
xlabel ("f [Hz]");
ylabel ("gain");
print (hfb, "gaindb_octave.odg", "-depsc");




