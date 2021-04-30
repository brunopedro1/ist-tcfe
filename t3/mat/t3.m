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

n = 10.20 

A = v_in/n


% Envelope Detector 

t=linspace(0, 0.1, 1000);

R_env = 10000  % envelope detector resistance 
C_env = 10e-6  % envelope detector capacitor

vS_env = A*cos(w*t)

vOhr = zeros(1, length(t));
vO_env = zeros(1, length(t));

tOFF = 1/w * atan(1/w/R_env/C_env);

vOnexp = A*cos(w*tOFF)*exp(-(t-tOFF)/R_env/C_env);

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

Vt = 0.026;
Is = 10e-14;
eta = 1
Rreg = 1e3

n_diodes = 18;
VON = 0.70;

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

hf=figure();
plot(t*1000, vOhr, t*1000, vO_env, t*1000, vO_reg);
title("Output voltage");
xlabel ("t[ms]");
legend("rectified","envelope", "regulator");
print (hf,"vO_env&vOhr&vO_reg.eps", "-depsc");


% Deviations (vO - 12) 

hfc = figure();
title("Deviations from DC voltage 12V")
plot (t*1000,vO_reg-12, ";vo-12 (t);");
xlabel ("t[ms]")
ylabel ("v_O [Volts]")
legend('Location','northeast');
print (hfc, "deviation.eps", "-depsc");





