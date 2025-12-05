%% Dylan Qiu Final Qube Motor Design ME351-2
% Fall 2025 
close all; clc;

% Motor Parameters
R = 8.4;                   % Resistance [Ohm]
L = 0.00116;               % Inductance [H]
J_m = 4.65e-6;             % Motor inertia [kg-m^2]
K_m = 0.042;               % Motor torque constant [N-m/A]
K_b = 0.042;               % Back EMF constant [V/(rad/s)]
B = 0;   

% Disc parameters
m_disc = 0.053;         
r_disc = 0.0248;      
J_d = 0.5 * m_disc * r_disc^2;  
J = J_m + J_d;  

num = K_m;  % Numerator
den = [L*J, (L*B + R*J), (R*B + K_m*K_b)];  % Denominator
sys = tf(num, den);

figure(1);
bode(sys);
r = roots(den);
title('Bode Plot of Motor Velocity Control')
[Gm,Pm,Wcg,Wcp] = margin(sys);
fprintf('\nGm = %.4f, Pm = %.4f, Wcg = %.4f, Wcp = %.4f \n', Gm,Pm,Wcg,Wcp);

for idx = 1:length(r)
    poles = r(idx);
    fprintf('Poles %d: pole = %.4f\n', idx, poles);
end

% damping and wn 

[wn,zeta] = damp(sys);
damp(sys);

%feedback(sys,1);
figure(2);
step(sys);

figure(3);
%controlSystemDesigner('rlocus',sys);

Kp = 7.55 / K_m; % found 7.55 from root loci plot breakaway pt
ess_step = 1 / (1+num*Kp); % steady state error for a type zero system is 1/ 1+K_vel and K_vel is simply the overall gain = product Kp and plant gain
fprintf('\nEss_step %d: \n', ess_step);
Kp_sys = Kp * sys;
bode(sys,Kp_sys);
%title('Bode Plot of Motor Velocity Control and Compensated System')
legend on;
%% Position Control
% Gp = K / s(ts+1)(ts+1)
tau_e = L / R;                  % Electrical Time Constant
tau_m = (R * J) / (K_m * K_b);  % Mechanical Time Constant
K_plant = 1/ K_b;               % Plant Gain

fprintf('Electrical Time Constant (tau_e): %.6f s\n', tau_e);
fprintf('Mechanical Time Constant (tau_m): %.6f s\n', tau_m);
fprintf('Plant Gain (K): %.4f\n', K_plant);

s = tf('s');
P = K_plant / (s * (tau_m*s + 1) * (tau_e*s + 1));

figure(4)

PM_target= 60; 
Phase_target = -180 + PM_target; 

[mag, phase, w] = bode(P, {0.1, 100}); 
phase = squeeze(phase); 
mag = squeeze(mag);
wc_new = interp1(phase, w, Phase_target);
mag_at_wc = interp1(w, mag, wc_new);
Kp = 1 / mag_at_wc;
Kp_P = Kp*P;
fprintf('Frequency for 60 deg PM: %.4f rad/s\n', wc_new);
fprintf('Required Proportional Gain (Kp): %.4f\n', Kp);

%bode(P, Kp * P)
K_c = 1;
t_r = 0.01; 
PI = K_c * (s + 1 / t_r) /s ;
PI_P = PI*P;
bode(P, Kp_P, PI_P);
legend on;
title('Bode Plot of Motor Position Control, P controller and PI controller')
%% Nyquist Plot for position control
figure(5)
title('Nyquist Plot for the plant')
nyquist1(P);

figure(6)
%PI_P_cl = feedback(PI_P,1);
nyquist1(PI_P);
title('Nyquist Plot position control with PI')

%% CONTROLLER DESIGN USING controlSystemDesigner MATLAB GUI 
fprintf('\n------CONTROLLER DESIGN USING controlSystemDesigner MATLAB GUI------\n')

%controlSystemDesigner(P);
%title('System designer of the plant ')

%% Lead control
close all; clc;
OS_max = 20;  % max over shoot percent
zeta_min = sqrt((log(OS_max/100))^2 / (pi^2 + (log(OS_max/100))^2));
pm = 100 * zeta_min; 

[~, ~, ~, Wcp_uncomp] = margin(P);
Wcp_comp = 3 * Wcp_uncomp; % new crossover frequency
[mag_wcp, phase_wcp] = bode(P, Wcp_comp);
fprintf('Desired compensated crossover freq: %.4f rad/s\n', Wcp_comp);

p_max = pm - (180 + phase_wcp) + 35; %35 degree margin for error
phi_max = p_max * pi/180;
alpha = (1 - sin(phi_max))/(1 + sin(phi_max));
z_lead = Wcp_comp * sqrt(alpha);
p_lead = z_lead / alpha;
fprintf('mag_wcp %.4f, phase_wcp %.4f\n', mag_wcp, phase_wcp);
Kc = 1 / (mag_wcp * sqrt(alpha));
fprintf('Compensator gain Kc: %.4f\n', Kc);

C_lead = Kc * (s + z_lead) / (s + p_lead);
P_Comp = P * C_lead;
fprintf('\nLead Compensator: Gc(s) = %.4f * (s + %.4f)/(s + %.4f)\n', Kc, z_lead, p_lead);
controlSystemDesigner(P_Comp);
P_comp_dist = feedback(P_Comp, C_lead);

dist_error = 1 - dcgain(P_comp_dist); 
fprintf('\nError due to unit disturbance = %.4f\n',dist_error);
%% Lead and/or Lag control

z_lag = 0.5;
p_lag = .04;
z_lead = 32;
p_lead = 32*5; 


P_comp_comp = P * Kc * (s+z_lead)*(s+z_lag) / ((s+p_lead)*(s+p_lag));
%controlSystemDesigner(P_comp_comp);
[~, ~, ~, Wcp_comp_comp] = margin(P_comp_comp);

fprintf('\nLead and lag Compensator: Gc(s) = %.4f * (s + %.4f)(s + %.4f)/(s + %.4f)(s + %.4f)\n', Kc, z_lead,z_lag, p_lead,p_lag);
fprintf('Lead + Lag compensated crossover freq= %.4f rad/s\n', Wcp_comp_comp);