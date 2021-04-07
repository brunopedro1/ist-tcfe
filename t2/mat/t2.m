% T2
clear all
format long 
pkg load control
pkg load symbolic
output_precision(10)

 %Units for the values: V, mA, kOhm, mS and uF
 %Values:  

fp = fopen('../mat/values.txt',"r");
values = dlmread(fp,'\n');
R1 = values(1)
R2 = values(2)
R3 = values(3)
R4 = values(4)
R5 = values(5)
R6 = values(6)
R7 = values(7)
Vs = values(8)
Kb = values(10)
Kd = values(11)
C = values(9)
null=0;
 
 
%R1 = 1.04001336091 
%R2 = 2.04372276851 
%R3 = 3.11359737601 
%R4 = 4.17085404861 
%R5 = 3.02859283303 
%R6 = 2.070545767 
%R7 = 1.01835949725 
%Vs = 5.20102702949 
%C = 1.00460501759 
%Kb = 7.19043597753 
%Kd = 8.06397385506 


%values = dlmread('data.txt');
%R1 = values(2,4) * (10^3);
%R2 = values(3,3) * (10^3);
%R3 = values(4,3) * (10^3);
%R4 = values(5,3) * (10^3);
%R5 = values(6,3) * (10^3);
%R6 = values(7,3) * (10^3);
%R7 = values(8,3) * (10^3);
%Vs = values(9,3);
%C = values(10,3) * (10^-6);
%Kb = values(11,3) * (10^-3);
%Kd = values(12,3) * (10^3);
%save("-ascii","../doc/tabelaVal.tex","R1", "R2", "R3", "R4", "R5", "R6", "R7", "Vs", "C", "Kb", "Kd");
%R= [R1, R2, R3, R4, R5, R6, R7];


%-------------------------------------------------------------------------------
% 1) Nodal analysis (t<0)
%-------------------------------------------------------------------------------
      
     % V1, V2, V3, V5, V6, V7, V8
A1 = [1,0,0,0,0,0,0; ... % ---------------------------------------node 1

      0,-(1/R1)-(1/R2)-(1/R3),(1/R2),(1/R3),0,0,0; ... % ---------node 2
      
      0,Kb+(1/R2),-(1/R2),-Kb,0,0,0; ...  % ----------------------node 3 
      
      0,-Kb,0,(1/R5)+Kb,-(1/R5),0,0; ... % -----------------------node 6
      
      0,0,0,0,0,-(1/R6)-(1/R7),(1/R7); ... % ---------------------node 7
      
      0,(1/R3),0,-(1/R3)-(1/R4)-(1/R5),(1/R5),(1/R7),-(1/R7); ... %supernode 5,8
      
      0,0,0,1,0,(Kd/R6),-1] % ------------------------------------additional equation: dependent voltage source
      
     
B1 = [Vs;-(Vs/R1);0;0;0;0;0]

C1 = A1\B1



V1 = C1(1,1)
V2 = C1(2,1)
V3 = C1(3,1)
V5 = C1(4,1)
V6 = C1(5,1)
V7 = C1(6,1)
V8 = C1(7,1)

Ib = Kb*(V2-V5)

IR1 = (V1 - V2)*(1/R1)
IR2 = (V2 - V3)*(1/R2)
IR3 = (V5 - V2)*(1/R3)
IR4 = V5*(1/R4)
IR5 = (V6 - V5)*(1/R5)
IR6 = -V7*(1/R6)
IR7 = (V7 - V8)*(1/R7)

%diary "../doc/nodal_tab.tex"
%diary on

%printf('v(1) = %.11f\n', V1);
%printf('v(2) = %.11f\n', V2);
%printf('v(3) = %.11f\n', V3);
%printf('v(5) = %.11f\n', V5);
%printf('v(6) = %.11f\n', V6);
%printf('v(7) = %.11f\n', V7);
%printf('v(8) = %.11f\n', V8);

%diary off

printf('op_TAB_nodal1\n');
printf('$V_1$ = %f\n$V_2$ = %f\n$V_3$ = %f\n$V_5$ = %f\n$V_6$ = %f\n$V_7$ = %f\n$V_8$ = %f\n$I_b$ = %f\n$R_1[i]$ = %f\n$R_2[i]$ = %f\n$R_3[i]$ = %f\n$R_4[i]$ = %f\n$R_5[i]$ = %f\n$R_6[i]$ = %f\n$R_7[i]$ = %f\n', V1, V2, V3, V5, V6, V7, V8, Ib, IR1,IR2,IR3,IR4,IR5,IR6,IR7);
printf('op_END_nodal1\n');


fp=fopen('../doc/ngspice_t21.tex',"w");
fprintf(fp,"v1 & %f \\\\ \\hline\nv2 & %f \\\\ \\hline\nv3 & %f \\\\ \\hline\nv4 & %f \\\\ \\hline\nv5 & %f \\\\ \\hline\nv6 & %f \\\\ \\hline\nv7 & %f \\\\ \\hline\nv8 & %f \\\\ \\hline",V1,V2,V3,null,V5,V6,V7,V8);
fclose(fp);

fp2 = fopen('../sim/ngspicevalues.txt',"w");
fprintf(fp2, "R1 1 2 %f\n",values(1)* (10^3));
fprintf(fp2, "R2 3 2 %f\n",values(2)* (10^3));
fprintf(fp2, "R3 2 5 %f\n",values(3)* (10^3));
fprintf(fp2, "R4 GND 5 %f\n",values(4)* (10^3));
fprintf(fp2, "R5 5 6 %f\n",values(5)* (10^3));
fprintf(fp2, "R6 GND 4 %f\n",values(6)* (10^3));
fprintf(fp2, "R7 7 8 %f\n",values(7)* (10^3));
fprintf(fp2, "C 6 8 %f\n",values(9)* (10^-6));
fprintf(fp2, "Vs 1 GND %f\n",values(8));
fprintf(fp2, "V3 4 7 DC 0 \n");
fprintf(fp2, "Hc 5 8 V3 %f\n",values(11)* (10^3));
fprintf(fp2, "Gb 6 3 (2,5) %f\n",values(10)* (10^-3));
fclose(fp2);

%-------------------------------------------------------------------------------
% 2) Nodal analysis (Vs = 0; Vx = V6-V8)
%-------------------------------------------------------------------------------
fp3 = fopen('../sim/ngspice_t22.txt',"w");
fprintf(fp3, "R1 GND 2 %f\n",values(1)* (10^3));
fprintf(fp3, "R2 3 2 %f\n",values(2)* (10^3));
fprintf(fp3, "R3 2 5 %f\n",values(3)* (10^3));
fprintf(fp3, "R4 GND 5 %f\n",values(4)* (10^3));
fprintf(fp3, "R5 5 6 %f\n",values(5)* (10^3));
fprintf(fp3, "R6 GND 4 %f\n",values(6)* (10^3));
fprintf(fp3, "R7 7 8 %f\n",values(7)* (10^3));
fprintf(fp3, "Vx 6 8 %f\n",V6-V8);
fprintf(fp3, "Vs 1 GND 0 \n");
fprintf(fp3, "V3 4 7 DC 0 \n");
fprintf(fp3, "Hc 5 8 V3 %f\n",values(11)* (10^3));
fprintf(fp3, "Gb 6 3 (2,5) %f\n",values(10)* (10^-3));

fclose(fp3);


Vs_2 = 0
V6 = C1(5,1)
V8 = C1(7,1)
Vx = V6-V8

    % V1, V2, V3, V5, V6, V7, V8, Ix
A2 = [1,0,0,0,0,0,0,0; ... %---------------------------------------node 1

     0,-(1/R1)-(1/R2)-(1/R3),(1/R2),(1/R3),0,0,0,0; ... % ---------node 2
     
     0,Kb+(1/R2),-(1/R2),-Kb,0,0,0,0; ...  % ----------------------node 3 
     
     0,-Kb,0,(1/R5)+Kb,-(1/R5),0,0,-1; ... % ----------------------node 6
     
     0,0,0,0,0,-(1/R6)-(1/R7),(1/R7),0; ... % ---------------------node 7
     
     0,(1/R3),0,-(1/R3)-(1/R4)-(1/R5),(1/R5),(1/R7),-(1/R7),1; ... %supernode 5,8
     
     0,0,0,1,0,(Kd/R6),-1,0; ... % --------------------------------additional equation: dependent voltage source
     
     0,0,0,0,1,0,-1,0] % ------------------------------------------additional equation: Vx = V6-V8
     
     
B2 = [Vs_2;-(Vs_2/R1);0;0;0;0;0;Vx]
     
C2 = A2\B2

V1_2 = C2(1,1)
V2_2 = C2(2,1)
V3_2 = C2(3,1)
V5_2 = C2(4,1)
V6_2 = C2(5,1)
V7_2 = C2(6,1)
V8_2 = C2(7,1)

Ix = C2(8,1)

Ib = Kb*(V2_2-V5_2)

IR1_2 = (V1_2 -V2_2)/R1
IR2_2 = (V3_2 -V2_2)/R2
IR3_2 = (V2_2 -V5_2)/R3
IR4_2 = V5_2/R4
IR5_2 = (V6_2 -V5_2)/R5
IR6_2 = -V7_2/R6
IR7_2 = (V7_2 -V8_2)/R7

Req = (Vx/Ix) % Equivalent resistance

tau = (Req*10^3)*(C*10^-6) % Time constant

printf('op_TAB_nodal2\n');
printf('$V_1$ = %f\n$V_2$ = %f\n$V_3$ = %f\n$V_5$ = %f\n$V_6$ = %f\n$V_7$ = %f\n$V_8$ = %f\n$I_b$ = %f\n$R_1[i]$ = %f\n$R_2[i]$ = %f\n$R_3[i]$ = %f\n$R_4[i]$ = %f\n$R_5[i]$ = %f\n$R_6[i]$ = %f\n$R_7[i]$ = %f\n$I_x$ = %f\n$R_{eq}$ = %f\n$tau$ = %f\n', V1_2, V2_2, V3_2, V5_2, V6_2, V7_2, V8_2, Ib, IR1_2,IR2_2,IR3_2,IR4_2,IR5_2,IR6_2,IR7_2,Ix, Req, tau);
printf('op_END_nodal2\n');

%-------------------------------------------------------------------------------
% 3) Natural solution, v_6n(t)
%-------------------------------------------------------------------------------
fp4 = fopen('../sim/ngspice_t23.txt',"w");
fprintf(fp4, "R1 1 2 %f\n",values(1)* (10^3));
fprintf(fp4, "R2 3 2 %f\n",values(2)* (10^3));
fprintf(fp4, "R3 2 5 %f\n",values(3)* (10^3));
fprintf(fp4, "R4 GND 5 %f\n",values(4)* (10^3));
fprintf(fp4, "R5 5 6 %f\n",values(5)* (10^3));
fprintf(fp4, "R6 GND 4 %f\n",values(6))* (10^3);
fprintf(fp4, "R7 7 8 %f\n",values(7)* (10^3));
fprintf(fp4, "C 6 8 %f\n",values(9)* (10^-6));
fprintf(fp4, "Vs 1 GND 0 \n");
fprintf(fp4, "V3 4 7 DC 0 \n");
fprintf(fp4, "Hc 5 8 V3 %f\n",values(11)* (10^3));
fprintf(fp4, "Gb 6 3 (2,5) %f\n",values(10)* (10^-3));

fclose(fp4);


t=0:0.000001:0.020;

v_6n = Vx*exp(-(t/abs(tau)));

hf = figure();
plot(t*10^3, v_6n, "r")
xlabel ("time [ms]");
ylabel ("V [V]");
grid on;
legend("v_6_n(t)");
print (hf, "natural_tab.eps", "-depsc");

%-------------------------------------------------------------------------------
% 4) Forced solution, v_6f(t)
%-------------------------------------------------------------------------------

f = 1000 % Hz
w= 2*pi*f
Zc = 1/(j*w*C)
Vs_4 = j

     % V1, V2, V3, V5, V6, V7, V8
A4 = [1,0,0,0,0,0,0; ... % ----------------------------------------------------node 1

      0,-(1/R1)-(1/R2)-(1/R3),(1/R2),(1/R3),0,0,0; ... % ----------------------node 2
      
      0,Kb+(1/R2),-(1/R2),-Kb,0,0,0; ...  % -----------------------------------node 3 
      
      0,-Kb,0,(1/R5)+Kb,-(1/R5)-(1/Zc),0,(1/Zc); ... % ------------------------node 6
      
      0,0,0,0,0,-(1/R6)-(1/R7),(1/R7); ... % ----------------------------------node 7
      
      0,(1/R3),0,-(1/R3)-(1/R4)-(1/R5),(1/Zc)+(1/R5),(1/R7),-(1/R7)-(1/Zc); ...%supernode 5,8
      
      0,0,0,1,0,(Kd/R6),-1] % -------------------------------------------------additional equation: dependent voltage source
      
     
B4 = [Vs_4;-(Vs_4/R1);0;0;0;0;0]

C4 = A4\B4

V1_f = C4(1,1)
V2_f = C4(2,1)
V3_f = C4(3,1)
V5_f = C4(4,1)
V6_f = C4(5,1)
V7_f = C4(6,1)
V8_f = C4(7,1)

magv1 = abs(V1_f)
magv2 = abs(V2_f)
magv3 = abs(V3_f)
magv5 = abs(V5_f)
magv6 = abs(V6_f)
magv7 = abs(V7_f)
magv8 = abs(V8_f)

printf('op_TAB_nodal4\n');
printf('$V_1$ = %f\n$V_2$ = %f\n$V_3$ = %f\n$V_5$ = %f\n$V_6$ = %f\n$V_7$ = %f\n$V_8$ = %f\n', magv1, magv2, magv3, magv5, magv6, magv7, magv8);
printf('op_END_nodal4\n');

phv1 = angle(V1_f)
phv2 = angle(V2_f)
phv3 = angle(V3_f)
phv5 = angle(V5_f)
phv6 = angle(V6_f)
phv7 = angle(V7_f)
phv8 = angle(V8_f)

printf('op_TAB_nodal4ph\n');
printf('$ \\phi_{V_1}$ = %f\n$\\phi_{V_2}$ = %f\n$\\phi_{V_3}$ = %f\n$\\phi_{V_5}$ = %f\n$\\phi_{V_6}$ = %f\n$\\phi_{V_7}$ = %f\n$\\phi_{V_8}$ = %f\n', phv1, phv2, phv3, phv5, phv6, phv7, phv8);
printf('op_END_nodal4ph\n');

t=0:0.000001:0.020;

v_6f(t>=0) = magv6*cos(w*t(t>=0)-phv6);
%v_6f(t>=0) = abs(V6_f)*sin(2*pi*f*t(t>=0));

hf2 = figure();
plot(t*10^3, v_6f, 'b')
xlabel ("time [ms]");
ylabel ("v_6_f [V]");
grid on;
legend("v_6_f(t)");
print (hf2, "force_tab.eps", "-depsc");


fp5 = fopen('../sim/ngspice_t24.txt',"w");
fprintf(fp5, "R1 1 2 %f\n",values(1)* (10^3));
fprintf(fp5, "R2 3 2 %f\n",values(2)* (10^3));
fprintf(fp5, "R3 2 5 %f\n",values(3)* (10^3));
fprintf(fp5, "R4 GND 5 %f\n",values(4)* (10^3));
fprintf(fp5, "R5 5 6 %f\n",values(5)* (10^3));
fprintf(fp5, "R6 GND 4 %f\n",values(6)* (10^3));
fprintf(fp5, "R7 7 8 %f\n",values(7)* (10^3));
fprintf(fp5, "C 6 8 %f\n",values(9)* (10^-6));
fprintf(fp5, "Vs 1 GND sin(0 1 1k) \n");
fprintf(fp5, "V3 4 7 DC 0 \n");
fprintf(fp5, "Hc 5 8 V3 %f\n",values(11)* (10^3));
fprintf(fp5, "Gb 6 3 (2,5) %f\n",values(10)* (10^-3));
fclose(fp5);

%-------------------------------------------------------------------------------
% 5) Total solution, v_6(t)
%-------------------------------------------------------------------------------

t=(-0.005:0.000001:0.020);

v_6n(t>=0) = Vx*exp(-(t(t>=0)/abs(tau)));

v_6f(t>=0) = magv6*cos(w*t(t>=0)-phv6);

v_6(t>=0) = v_6n(t>=0) + v_6f(t>=0);
v_6(t<0) = V6;

v_s(t>=0) = sin(2*pi*f*t(t>=0));
v_s(t<0) = Vs;


hf3 = figure();
plot(t*1000, v_6, t*1000, v_s)
xlabel ("time [ms]");
ylabel ("voltage [V]");
grid on;
legend("v_6(t)", "v_s(t)");
print (hf3, "total_tab.eps", "-depsc");

%-------------------------------------------------------------------------------
% 6) Frequency Response
%-------------------------------------------------------------------------------

fp6 = fopen('../sim/ngspice5.txt',"w");
fprintf(fp6, "R1 1 2 %f\n",values(1)* (10^3));
fprintf(fp6, "R2 3 2 %f\n",values(2)* (10^3));
fprintf(fp6, "R3 2 5 %f\n",values(3)* (10^3));
fprintf(fp6, "R4 GND 5 %f\n",values(4)* (10^3));
fprintf(fp6, "R5 5 6 %f\n",values(5)* (10^3));
fprintf(fp6, "R6 GND 4 %f\n",values(6)* (10^3));
fprintf(fp6, "R7 7 8 %f\n",values(7)* (10^3));
fprintf(fp6, "C 6 8 %f\n",values(9)* (10^-6));
fprintf(fp6, "Vs 1 GND 1 ac 1 sin(0 1 1k) \n");
fprintf(fp6, "V3 4 7 DC 0 \n");
fprintf(fp6, "Hc 5 8 V3 %f\n",values(11)* (10^3));
fprintf(fp6, "Gb 6 3 (2,5) %f\n",values(10)* (10^-3));
fclose(fp6);


fp = fopen('../mat/values.txt',"r");
values = dlmread(fp,'\n');
R16 = values(1)* (10^3)
R26 = values(2)* (10^3)
R36 = values(3)* (10^3)
R46 = values(4)* (10^3)
R56 = values(5)* (10^3)
R66 = values(6)* (10^3)
R76 = values(7)* (10^3)
Vs6 = values(8)
Kb6 = values(10)* (10^-3)
Kd6 = values(11)* (10^3)
C6 = values(9)* (10^-6)
null=0;


f = logspace(-1,6,200);


Vs = 1;

%Zc1 = 1/(j*w*(C*(10^-6)))


for i= 1:length(f)

w = 2*pi*f(i);
Zc1 = 1/(j*w*C6);

U = [1/(R16),-(1/(R16))-(1/(R26))-(1/(R36)),1/(R26),1/(R36),0,0,0; 0,((Kb6)+1/(R26)),-1/(R26),-(Kb6),0,0,0;0,(Kb6),0,(-(Kb6)-(1/(R56))),(1/(R56)+1/Zc1),0,-(1/Zc1); 1,0,0,0,0,0,0; 0,0,0,1,0,((Kd6)/(R66)),-1; 0,0,0,0,0,(1/(R66) +1/(R76)), -(1/(R76)); 0,1/(R36),0,-(1/(R36))-(1/(R46))-(1/(R56)),(1/(R56) +1/Zc1),1/(R76), (-(1/(R76))-(1/Zc1))];

Y = inv(U);
W = [0;0;0;Vs;0;0;0];	

K = Y*W;

v6(i)= K(5);
vs(i)= K(1);
vc(i)= K(5)-K(7);


endfor


for i = 1:length(f)

v6amp(i)= 20*log10(abs(v6(i)));
v6phs(i)= (180/pi)*angle(v6(i));
vsamp(i)= 20*log10(abs(vs(i)));
vsphs(i)= (180/pi)*angle(vs(i));
vcamp(i)= 20*log10(abs(vc(i)));
vcphs(i)= (180/pi)*angle(vc(i));



endfor



f6 = figure();

semilogx(f,v6amp,"r");
hold on

semilogx(f,vcamp,"g");
hold on

semilogx(f,vsamp,"c");
hold on

 

xlabel('f[Hz]');

ylabel('V[dB]');

title("Frequency Response- Amplitude");

legend('v6amp','vcamp','vsamp','Location','Northeast');

print(f6,"FrequencyResponseAmplitude.eps","-depsc");






f7 = figure();

semilogx(f,v6phs,"r");
hold on

semilogx(f,vcphs,"g");
hold on

semilogx(f,vsphs,"c");
hold on


xlabel('f[Hz]');

ylabel('Angle[degrees]');

title("Frequency Response- Phase");

legend('v6phs','vcphs','vsphs','Location','Northeast');

print(f7,"FrequencyResponsePhase.eps","-depsc");
