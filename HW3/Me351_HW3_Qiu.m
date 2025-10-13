% ME351 HW3 Matlab Code
% Dylan Qiu, ME '27

close all; clc;

%% PROBLEM 1

fprintf('PROBLEM 1\n');
% Transfer function: G1(s) = K / (s(s^2+2s+2))
K = 0.1;
s = tf('s');
G1 = K / (s * (s^2 + 2 * s + 2));

%controlSystemDesigner('rlocus', G1)

%% PROBLEM 2

fprintf('PROBLEM 2\n');
% Transfer function: G1(s) = K / (s(s+2)(s^2+2s+2))
K = 1;
s = tf('s');
G2 = K / (s * (s+2) * (s^2 + 2 * s + 2));

%controlSystemDesigner('rlocus', G2)

%% PROBLEM 3

fprintf('PROBLEM 3\n');
% Transfer function: G1(s) = K(s+2) / (s(s^2+2s+2))
K = 0.01;
s = tf('s');
G3 = (K * (s+2)) / (s * (s^2 + 2 * s + 2));

%controlSystemDesigner('rlocus', G3)

%% PROBLEM 4

fprintf('PROBLEM 4a\n');
% Transfer function: G1(s) = K / (s^2)
K = 0.01;
J = 1 ; %kg-m^2
s = tf('s');
G4a = (K)/ (J*s^2);

%controlSystemDesigner('rlocus', G4)


fprintf('PROBLEM 4g\n');
% Transfer function: G1(s) = s^2 + K(s+1)
K = 0.01;
J = 1 ; %kg-m^2
s = tf('s');
G4g = K * (s + 1) / s^2 ;
%controlSystemDesigner('rlocus', G4g)

z_lead = 1;
p_values = [20, 9, 3];

%figure('Problem 4g\n');

for idx = 1:length(p_values)
    p_lead = p_values(idx);
    
    fprintf('\n(g.%d) p = %d:\n', idx, p_lead);
    
    % Transfer function
    sg = tf('s');
    G_lead = (sg + 1) / (sg^3 + p_lead * sg^2);
    G_lead;
    
    %subplot(1, 3, idx);
    %rlocus(G_lead);
    %grid on;
    %title(sprintf('p = %d', p_lead));
    %xlabel('Real Axis');
    %ylabel('Imaginary Axis');
    %sgrid;
    %axis equal;
    %xlim([-p_lead-2, 2]);
end

%% PROBLEM 5
fprintf('Problem 5:');

% Given parameters
R = 8.4;                    % Resistance (Ohms)
L = 0.00116;                % Inductance (H)
Jm = 4.65e-6;               % Motor inertia (kg-m^2)
Km = 0.042;                 % Motor torque constant (N-m/A)
Kb = 0.042;                 % Back-EMF constant (V/(rad/s))
B = 0.02;                   % Viscous damping (N-m-s)

% Disc parameters
m_disc = 0.053;             % Disc mass (kg)
r_disc = 0.0248;            % Disc radius (m)
Jd = 0.5 * m_disc * r_disc^2;  % Disc inertia (kg-m^2)

% Total inertia
J_total = Jm + Jd;

% Transfer function coefficients
a2 = L * J_total;
a1 = R * J_total + L * B;
a0 = R * B + Km * Kb;

% Open-loop transfer function
num5 = Km;
den5 = [a2, a1, a0];
G5 = tf(num5, den5);

fprintf('Transfer Function G(s) = omega(s)/V(s):\n');
disp(G5);

% Open-loop poles
poles5 = pole(G5);
fprintf('Open-Loop Poles:\n');
disp(poles5);

% Design specifications
tr_spec = 0.5e-3;           % Rise time <= 0.5 ms
Mp_spec = 0.05;             % Overshoot <= 5%

% Calculate requirements
zeta_min = sqrt((log(Mp_spec))^2 / (pi^2 + (log(Mp_spec))^2));
wn_min = 1.8 / (tr_spec);
zeta_min;
wn_min;

%% Part (a)

controlSystemDesigner('rlocus', G5);

%% Part (f): Minimum controller gain
fprintf('\n--- Part (f): Minimum Controller Gain ---\n');
min_gain_ratio = ((1 - 0.02) / 0.02) / dcgain(G5);
min_gain_ratio;






