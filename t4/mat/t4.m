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

%op
RB=1/(1/RB1+1/RB2)
VEQ=RB2/(RB1+RB2)*VCC
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1)
IC1=BFN*IB1
IE1=(1+BFN)*IB1
VE1=RE1*IE1
VO1=VCC-RC1*IC1
VCE=VO1-VE1


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

RE1=100
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)))
ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) ) 
ZO1 = 1/(1/ZX+1/RC1)

RE1=0
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZO1 = 1/(1/ro1+1/RC1)

%ouput stage
BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2

AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)


%total
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)




f         = logspace(1, 6+2, 10*(log10(100e6/10)));
w         = f*2*pi;
GainFreq  = w;


for c=1:size(w,2)
  
  ZCI = 1/(w(c)*j*CI);
  ZE1 = RE1;
  ZCB = 1/(w(c)*j*CB);
  ZE1CB = 1/(1/ZE1 + 1/ZCB);
  ZCO = 1/(w(c)*j*CO);
  ZL = RL
  ZS = RS
  
  AG = [  ZS+ZCI+ZB , -ZB               , 0   , 0               , 0           , 0     , 0             ; ...
        -ZB       , ZB + zpi1 + ZE1CB , 0   , -ZE1CB          , 0           , 0     , 0             ; ...
        0         , zpi1*gm1          , 1   , 0               , 0           , 0     , 0             ; ...
        0         , -ZE1CB            , -zo1, ZE1CB + zo1 + ZC, -ZC         , 0     , 0             ; ...
        0         , 0                 , 0   , -ZC             , zpi2+ZC     , Zo2E2 , -Zo2E2        ; ...
        0         , 0                 , 0   , 0               , -1-zpi2*gm2 , 1     , 0             ; ...
        0         , 0                 , 0   , 0               , 0           , -Zo2E2, Zo2E2+ZCO+ZL ];
        
  BG = [  Vin       ; 0                 ; 0   ; 0               ; 0           ; 0     ; 0     ];

  PG = AG\BG;

  GainFreq(c) = PG(7)*ZL/Vin;
  
endfor

f_H=f;
f_L= (1/(3*CS)+ 1/((ZO+RL)*Co) +1/((ZI+RS)*Cb))/(2*pi);
band=f_H - f_L

freq=logspace(1,8);

for i=1:50
wfreq=2*pi*freq(i)
ZCS   = 1/(wfreq*j*CS);
ZCB   = 1/(wfreq*j*Cb);
ZEB = 1/(1/RE1 + 1/ZCB);
zpi2  = 1/gpi2;
zo2   = 1/go2;
ZE2 = 1/(1/zo2 + 1/RE2);
ZCO   = 1/(wfreq*j*Co);

 HOT = [  RS+ZCS+RB , -RB, 0 , 0 , 0 , 0 , 0; ...
 	0 , -ZEB , -ro1, ZEB + ro1 + RC1, -RC1, 0, 0 ; ...
        0 , rpi1*gm1 , 1  , 0 , 0, 0 , 0  ; ...
        0 , 0 , 0  , -RC1, zpi2+RC1, ZE2 , -ZE2; ...
        0 , 0  , 0  , 0  , -1-zpi2*gm2 , 1  , 0 ; ...
        0, 0 , 0  , 0 , 0 , -ZE2, ZE2+ZCO+RL ; ...
        -RB  , RB + rpi1 + ZEB , 0 , -ZEB  , 0 , 0 , 0 ];
        
DOG = [Vinput; 0; 0 ; 0 ; 0; 0; 0];

Gain= HOT\DOG;

GaindB(i) = 20*log10(abs(Gain(7)*RL/Vinput))
endfor






