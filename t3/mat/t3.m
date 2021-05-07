close all
clear all
format long

%-------------------------------------------------------------------------------
% AC/DC CONVERTER
%-------------------------------------------------------------------------------

f = 50
w=2*pi*f
v_in = 230


% Transformer n:1

n = 10.7785890358

A = v_in/n


% Envelope Detector 

t=linspace(0, 10/f, 1000);

R_env = 150000  % envelope detector resistance 
C_env = 100e-6  % envelope detector capacitor

vS_env = A*cos(w*t)

vOhr = zeros(1, length(t));
vO_env = zeros(1, length(t));

tOFF = 1/w * atan(1/w/R_env/C_env);

vOnexp = A*cos(w*tOFF)*exp(-(t-tOFF)/R_env/C_env);

%Half wave rectifier
%for i=1:length(t)
%	  if (vS_env(i) > 0)
%	    vOhr(i) = vS_env(i);
%	  else
%	    vOhr(i) = 0;
%	  endif
%endfor

%Full wave rectifier
for i=1:length(t)
	  vOhr(i) = abs(vS_env(i));
endfor


for i=1:length(t)
  if t(i) < tOFF
    vO_env(i) = vOhr(i);
  elseif vOnexp(i) > vOhr(i) % t = tON
    vO_env(i) = vOnexp(i);
  else 
    tOFF = tOFF + 1/f/2;
    vOnexp = A*abs(cos(w*tOFF))*exp(-(t-tOFF)/R_env/C_env);
    vO_env(i) = vOhr(i);
  endif
endfor

ripple_env = max(vO_env) - min(vO_env)



% Voltage Regulator

vO_reg = zeros(1, length(t));
VO_reg = 0;
vo_reg = zeros(1, length(t));

Vt = 0.026;
Is = 1e-14;
eta = 1
Rreg = 5000

n_diodes = 18;
VON = 12/n_diodes;

Rd = Vt*eta/(Is*exp(VON/(eta*Vt)))

if mean(vO_env) >= VON*n_diodes
  VO_reg = VON*n_diodes;
else
  VO_reg = mean(vO_env);
endif


for i = 1:length(t)
  if vO_env(i) >= n_diodes*VON
    vo_reg(i) = n_diodes*Rd/(n_diodes*Rd+Rreg) * (vO_env(i)-mean(vO_env));
  else
    vo_reg(i) = vO_env(i)-mean(vO_env);
  endif
endfor

vO_reg = VO_reg + vo_reg;

ripple_reg = max(vO_reg)-min(vO_reg) 

%hf=figure();
%plot(t*1000, vOhr, t*1000, vO_env, t*1000, vO_reg);
%title("Output voltage");
%xlabel ("t[ms]");
%legend("rectified","envelope", "regulator");
%print (hf,"vO_env&vOhr&vO_reg.eps", "-depsc");

%hfr=figure();
%plot(t*1000, vOhr);
%xlabel ("t [ms]");
%ylabel ("V [Volts]");
%legend("V_{rectified}");
%print (hfr,"vOhr.eps", "-depsc");

hfenv=figure();
plot(t*1000, vO_env);
xlabel ("t [ms]");
ylabel ("V [Volts]");
legend("V_{envelope}");
print (hfenv,"vO_env.eps", "-depsc");

hfreg=figure();
plot(t*1000, vO_reg);
xlabel ("t [ms]");
ylabel ("V [Volts]");
legend("V_{regulator}");
print (hfreg,"vO_reg.eps", "-depsc");

% Deviations (vO - 12) 

hfc = figure();
title("Deviations from DC voltage 12V")
plot (t*1000,vO_reg-12);
xlabel ("t[ms]")
ylabel ("V [Volts]")
legend('V_{O} - 12 ');
print (hfc, "deviation.eps", "-depsc");


cost = R_env/1000 + Rreg/1000 + C_env*1e6 + n_diodes*0.1+0.4

merit = 1/(cost*(ripple_reg + abs(mean(vO_reg) - 12) + 1e-6))

printf('op_TAB\n');
printf('n = %f\n$R_{envelope}$  [kOhm] = %f\n$C_{envelope}$  [uF] = %f\n$R_{regulator}$  [kOhm] = %f\n', n, R_env/1000, C_env*1e6, Rreg/1000);
printf('op_END\n');

printf('op2_TAB\n');
printf('$Avarage(v_O)$  [V] = %f\n$Ripple(v_O)$  [V] = %f\n$Cost$  [MU] = %f\n$Merit$ = %f\n', mean(vO_reg), ripple_reg , cost, merit);
printf('op2_END\n');
