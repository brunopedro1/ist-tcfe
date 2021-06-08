close all
clear all
format long

R1 = 1000
R2 = 1000
R3 = 100000
R4 = 1000
C1 = 220e-9
C2 = 220e-9

printf("UsedValues_TAB \n"); 
printf("$R_{1}$ = %e Ohm \n", R1);
printf("$R_{2}$ = %e Ohm \n", R2);
printf("$R_{3}$ = %e Ohm \n", R3);
printf("$R_{4}$ = %e Ohm \n", R4);
printf("$C_{1}$ = %e F \n", C1);
printf("$C_{2}$ = %e F \n", C2);
printf("UsedValues_END \n \n");

%central frequency
wL = 1/(R1*C1)
wH = 1/(R2*C2)
wO = sqrt(wL*wH)
f = wO/(2*pi)

%gain
gain = abs((R1*C1*wO*j)/(1+R1*C1*wO*j)*(1+R3/R4)*(1/(1+R2*C2*wO*j)))
gaindb = 20*log10(abs(gain))

%input impedance
Zin = abs(R1 + 1/(j*wO*C1))

%output impedance
Zout = abs(1/(j*wO*C2+1/R2))

printf("teoresults_TAB \n"); 
printf("Central frequency = %e Hz \n", f);
printf("Gain = %e dB\n", gaindb);
printf("$Z_{input}$ = %e Ohm \n", Zin);
printf("$Z_{output}$ = %e Ohm \n", Zout);
printf("teoresults_END \n \n");

%frequency response
freq = logspace(1,8,100);

rp = (R1*C1*2*pi*freq*j)./(1+R1*C1*2*pi*freq*j)*(1+R3/R4).*(1./(1+R2*C2*2*pi*freq*j));

f1 = figure();
semilogx(freq,20*log10(abs(rp)));
xlabel("Frequency [Hz]");
ylabel("Gain [dB]");
title("Gain");
print(f1, "teo_gain.eps", "-depsc");

f2 = figure();
semilogx(freq,180*arg(rp)/pi);
xlabel("Frequency [Hz]");
ylabel("Phase [Deg]");
title("Phase");
print(f2, "teo_phase.eps", "-depsc");
