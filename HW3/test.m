%% Problem 5 

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
a0 = R * B + Km * Kb

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
zeta_min
wn_min

%% Part (a)

controlSystemDesigner('rlocus', G5);

%% Part (f): Minimum controller gain
fprintf('\n--- Part (f): Minimum Controller Gain ---\n');
min_gain_ratio = ((1 - 0.02) / 0.02) / dcgain(G5);
min_gain_ratio




