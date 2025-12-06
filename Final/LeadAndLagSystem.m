%% Position Control Lead and Lag Compensator
% Gp = K / s(ts+1)(ts+1)

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

tau_e = L / R;                  % Electrical Time Constant
tau_m = (R * J) / (K_m * K_b);  % Mechanical Time Constant
K_plant = 1/ K_b;               % Plant Gain

%fprintf('Electrical Time Constant (tau_e): %.6f s\n', tau_e);
%fprintf('Mechanical Time Constant (tau_m): %.6f s\n', tau_m);
fprintf('Plant Gain (K): %.4f\n', K_plant);

s = tf('s');
P = K_plant / (s* (tau_m*s + 1) * (tau_e*s + 1));

OS_max = 20;  % max over shoot percent
zeta_min = sqrt((log(OS_max/100))^2 / (pi^2 + (log(OS_max/100))^2));
pm = 100 * zeta_min; 

[~, ~, ~, Wcp_uncomp] = margin(P);
Wcp_comp = 3 * Wcp_uncomp; % new crossover frequency
[mag_wcp, phase_wcp] = bode(P, Wcp_comp);
fprintf('Desired compensated crossover freq: %.4f rad/s\n', Wcp_comp);

p_max = pm - (180 + phase_wcp) + 30; % 40 degree margin for error
phi_max = p_max * pi/180; % convert to radians
alpha = (1 - sin(phi_max))/(1 + sin(phi_max));
z_lead = Wcp_comp * sqrt(alpha);
p_lead = z_lead / alpha;
fprintf('mag_wcp: %.4f, phase_wcp: %.4f\n', mag_wcp, phase_wcp);
Kc = 1 / (mag_wcp * sqrt(alpha));
fprintf('Compensator gain Kc: %.4f\n', Kc);

z_lag = 4;
p_lag = z_lag/75; % for lag 
C_lead = (s + z_lead) / (s + p_lead);
C_lag  = (s + z_lag)  / (s + p_lag);
C = Kc * C_lead * C_lag;
P_comp_comp = P * C;
controlSystemDesigner(P); %Kc = 104.7  zeroes: (s+35) (s+0.06) poles:   (s+165) (s+0.028)


fprintf('\nLead and lag Compensator: Gc(s) = %.4f * (s + %.4f)(s + %.4f)/(s + %.4f)(s + %.4f)\n', Kc, z_lead,z_lag, p_lead,p_lag);
P_dist = feedback(P_comp_comp,C);
P_step = feedback(P_comp_comp,1);

info = stepinfo(P_step);

% Calculate Steady State Error for Disturbance (1V Step)
dist_error = dcgain(P_dist); 

% Reference Steady State Error (Step Input)
ref_error = 1 - dcgain(P_step);

fprintf('\nDESIGN RESULTS:\n');
fprintf('1. Rise Time (Target <= 0.1s)=      %.4f s\n', info.RiseTime);
fprintf('2. Overshoot (Target < 20%%)=     %.2f %%\n', info.Overshoot);
fprintf('3. Step input SS Error (Target ~0)=     %.4f\n', ref_error);
fprintf('4. Disturbance Error (Target < 0.05)= %.4f rad\n', dist_error);

