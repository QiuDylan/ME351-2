% DC Motor To Sinusoidal Input Analysis
% Authors: Dylan Qiu
% Date: 10/22/2025

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
den = [L*J, (L*B + R*J), (R*B + K_m*K_b)];  % Denominator coefficients
sys = tf(num, den);

grid on;
figure(1);
step(sys);

grid on;
figure(2);
bode(sys);