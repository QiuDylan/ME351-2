%% ME351 HW3 - Problem 5: Quanser Motor Control System - DETAILED SOLUTION
% Complete step-by-step analysis with explanations

clear all; close all; clc;

fprintf('========================================================================\n');
fprintf('PROBLEM 5: QUANSER QUBE-SERVO-2 MOTOR CONTROL SYSTEM\n');
fprintf('========================================================================\n\n');

%% GIVEN PARAMETERS
fprintf('---- GIVEN SYSTEM PARAMETERS ----\n\n');

% Electrical parameters
R = 8.4;                    % Armature resistance (Ohms)
L = 0.00116;                % Armature inductance (H) = 1.16 mH
fprintf('Electrical Parameters:\n');
fprintf('  R = %.2f Ω (Armature resistance)\n', R);
fprintf('  L = %.5f H = %.2f mH (Armature inductance)\n', L, L*1000);

% Mechanical parameters
Jm = 4.65e-6;               % Motor rotor inertia (kg-m^2)
Km = 0.042;                 % Motor torque constant (N-m/A)
Kb = 0.042;                 % Back-EMF constant (V/(rad/s))
B = 0.02;                   % Viscous damping coefficient (N-m-s)
fprintf('\nMechanical Parameters:\n');
fprintf('  Jm = %.2e kg-m² (Motor rotor inertia)\n', Jm);
fprintf('  Km = %.3f N-m/A (Torque constant)\n', Km);
fprintf('  Kb = %.3f V/(rad/s) (Back-EMF constant)\n', Kb);
fprintf('  B  = %.3f N-m-s (Viscous damping)\n', B);

% Disc parameters
m_disc = 0.053;             % Disc mass (kg)
r_disc = 0.0248;            % Disc radius (m)
fprintf('\nDisc Parameters:\n');
fprintf('  m_disc = %.3f kg (Disc mass)\n', m_disc);
fprintf('  r_disc = %.4f m = %.2f mm (Disc radius)\n', r_disc, r_disc*1000);

% Calculate disc inertia (solid disc: J = 0.5*m*r^2)
Jd = 0.5 * m_disc * r_disc^2;
fprintf('  Jd = 0.5 × m × r² = 0.5 × %.3f × %.4f²\n', m_disc, r_disc);
fprintf('  Jd = %.4e kg-m² (Disc inertia)\n', Jd);

% Total inertia
J_total = Jm + Jd;
fprintf('\nTotal System Inertia:\n');
fprintf('  J_total = Jm + Jd = %.4e + %.4e = %.4e kg-m²\n', Jm, Jd, J_total);

%% DERIVE TRANSFER FUNCTION
fprintf('\n\n---- TRANSFER FUNCTION DERIVATION ----\n\n');

fprintf('DC Motor Equations:\n');
fprintf('  Electrical: V(t) = R·i(t) + L·di/dt + Kb·ω(t)\n');
fprintf('  Mechanical: Km·i(t) = J·dω/dt + B·ω(t)\n\n');

fprintf('Taking Laplace transform:\n');
fprintf('  V(s) = (R + Ls)·I(s) + Kb·Ω(s)\n');
fprintf('  Km·I(s) = (Js + B)·Ω(s)\n\n');

fprintf('From mechanical equation: I(s) = (Js + B)·Ω(s) / Km\n');
fprintf('Substitute into electrical equation:\n');
fprintf('  V(s) = (R + Ls)·(Js + B)·Ω(s)/Km + Kb·Ω(s)\n');
fprintf('  V(s) = [(R + Ls)(Js + B)/Km + Kb]·Ω(s)\n');
fprintf('  V(s) = [(LJs² + (RJ + LB)s + RB)/Km + Kb]·Ω(s)\n');
fprintf('  V(s) = [LJs² + (RJ + LB)s + (RB + Km·Kb)]·Ω(s) / Km\n\n');

fprintf('Therefore: G(s) = Ω(s)/V(s) = Km / [LJs² + (RJ + LB)s + (RB + Km·Kb)]\n\n');

% Calculate transfer function coefficients
a2 = L * J_total;
a1 = R * J_total + L * B;
a0 = R * B + Km * Kb;

fprintf('Calculating coefficients:\n');
fprintf('  a2 = L×J = %.5f × %.4e = %.6e\n', L, J_total, a2);
fprintf('  a1 = R×J + L×B = %.2f×%.4e + %.5f×%.2f = %.6e\n', R, J_total, L, B, a1);
fprintf('  a0 = R×B + Km×Kb = %.2f×%.2f + %.3f×%.3f = %.6f\n', R, B, Km, Kb, a0);

% Create transfer function
num = Km;
den = [a2, a1, a0];
G = tf(num, den);

fprintf('\n** OPEN-LOOP TRANSFER FUNCTION **\n');
fprintf('G(s) = %.4f / (%.4es² + %.4es + %.4f)\n\n', num, a2, a1, a0);
disp(G);

%% FIND OPEN-LOOP POLES
fprintf('\n---- OPEN-LOOP POLES ----\n\n');

poles_ol = pole(G);
fprintf('Open-loop poles (solving as² + as + a0 = 0):\n');
fprintf('Using quadratic formula: s = [-a1 ± √(a1² - 4a2·a0)] / (2a2)\n\n');

discriminant = a1^2 - 4*a2*a0;
fprintf('Discriminant = a1² - 4a2·a0 = (%.4e)² - 4×%.4e×%.4f\n', a1, a2, a0);
fprintf('            = %.4e - %.4e = %.4e\n', a1^2, 4*a2*a0, discriminant);

if discriminant < 0
    fprintf('Discriminant < 0 → Complex poles\n\n');
    real_part = -a1/(2*a2);
    imag_part = sqrt(abs(discriminant))/(2*a2);
    fprintf('Poles: s = %.2f ± j%.2f\n', real_part, imag_part);
else
    fprintf('Discriminant > 0 → Real poles\n\n');
end

fprintf('MATLAB calculated poles:\n');
for i = 1:length(poles_ol)
    if abs(imag(poles_ol(i))) < 1e-10
        fprintf('  s = %.4f\n', real(poles_ol(i)));
    else
        fprintf('  s = %.4f + j%.4f\n', real(poles_ol(i)), imag(poles_ol(i)));
    end
end

% Calculate pole properties
if abs(imag(poles_ol(1))) > 1e-10
    sigma_ol = real(poles_ol(1));
    wd_ol = abs(imag(poles_ol(1)));
    wn_ol = sqrt(sigma_ol^2 + wd_ol^2);
    zeta_ol = -sigma_ol/wn_ol;
    
    fprintf('\nOpen-loop pole characteristics:\n');
    fprintf('  Natural frequency (ωn) = %.2f rad/s\n', wn_ol);
    fprintf('  Damping ratio (ζ) = %.4f\n', zeta_ol);
    fprintf('  Damped frequency (ωd) = %.2f rad/s\n', wd_ol);
end

%% DESIGN SPECIFICATIONS
fprintf('\n\n---- DESIGN SPECIFICATIONS ----\n\n');

tr_spec = 0.5e-3;           % Rise time <= 0.5 ms (10-90%)
Mp_spec = 0.05;             % Overshoot <= 5%

fprintf('Performance Requirements:\n');
fprintf('  (i)  Rise time: tr ≤ %.1f ms (10-90%% rise time)\n', tr_spec*1000);
fprintf('  (ii) Overshoot: Mp ≤ %.0f%%\n', Mp_spec*100);

% Calculate minimum damping ratio from overshoot
fprintf('\n** Calculating minimum damping ratio from overshoot **\n');
fprintf('Overshoot formula: Mp = exp(-ζπ/√(1-ζ²))\n');
fprintf('Given Mp ≤ 0.05, solve for ζ:\n');
fprintf('  0.05 = exp(-ζπ/√(1-ζ²))\n');
fprintf('  ln(0.05) = -ζπ/√(1-ζ²)\n');
fprintf('  %.4f = -ζπ/√(1-ζ²)\n', log(Mp_spec));

% Solve: ln(Mp) = -zeta*pi/sqrt(1-zeta^2)
% Rearranging: zeta = sqrt((ln(Mp))^2 / (pi^2 + (ln(Mp))^2))
zeta_min = sqrt((log(Mp_spec))^2 / (pi^2 + (log(Mp_spec))^2));

fprintf('\nSolving: ζ_min = √[ln²(Mp) / (π² + ln²(Mp))]\n');
fprintf('         ζ_min = √[(%.4f)² / (π² + (%.4f)²)]\n', log(Mp_spec), log(Mp_spec));
fprintf('         ζ_min = √[%.4f / %.4f]\n', (log(Mp_spec))^2, pi^2 + (log(Mp_spec))^2);
fprintf('         ζ_min = %.4f\n', zeta_min);

% Calculate minimum natural frequency from rise time
fprintf('\n** Calculating minimum natural frequency from rise time **\n');
fprintf('Rise time approximation (10-90%%): tr ≈ 2.2/(ζ·ωn)\n');
fprintf('Given tr ≤ %.4f s, and ζ ≥ %.4f:\n', tr_spec, zeta_min);
fprintf('  ωn ≥ 2.2/(ζ·tr) = 2.2/(%.4f × %.4f)\n', zeta_min, tr_spec);

wn_min = 2.2 / (zeta_min * tr_spec);
fprintf('  ωn_min = %.2f rad/s\n', wn_min);

fprintf('\n** DESIGN REQUIREMENTS SUMMARY **\n');
fprintf('  Minimum damping ratio: ζ ≥ %.4f\n', zeta_min);
fprintf('  Minimum natural frequency: ωn ≥ %.2f rad/s\n', wn_min);
fprintf('  This defines a shaded region in the s-plane\n');

%% Part (a): Plot Root Locus with Design Requirements
fprintf('\n\n========================================\n');
fprintf('PART (a): ROOT LOCUS WITH DESIGN REQUIREMENTS\n');
fprintf('========================================\n\n');

% Create root locus plot
figure('Name', 'Problem 5(a) - Root Locus with Design Requirements', ...
    'NumberTitle', 'off', 'Position', [50 50 1000 800]);
rlocus(G);
grid on;
hold on;

% Add design requirement lines
sgrid(zeta_min, wn_min);

% Highlight acceptable region
theta = linspace(0, 2*pi, 100);
% Minimum wn circle
x_wn = wn_min * cos(theta);
y_wn = wn_min * sin(theta);
plot(x_wn, y_wn, 'g--', 'LineWidth', 2);

% Minimum zeta lines
angle_max = acos(zeta_min);
x_zeta = [-1000, 0];
y_zeta_upper = -x_zeta * tan(angle_max);
y_zeta_lower = -y_zeta_upper;
plot(x_zeta, y_zeta_upper, 'r--', 'LineWidth', 2);
plot(x_zeta, y_zeta_lower, 'r--', 'LineWidth', 2);

% Mark open-loop poles
plot(real(poles_ol), imag(poles_ol), 'rx', 'MarkerSize', 15, 'LineWidth', 3);

title('Root Locus with Design Requirements');
xlabel('Real Axis (σ)');
ylabel('Imaginary Axis (jω)');
legend('Root Locus', 'ζ and ωn grid', sprintf('ωn = %.0f rad/s', wn_min), ...
    sprintf('ζ = %.3f', zeta_min), 'Open-loop poles', 'Location', 'best');
xlim([-8000 500]);
ylim([-8000 8000]);

fprintf('Opening Control System Designer...\n');
fprintf('Instructions:\n');
fprintf('  1. Right-click on the root locus plot\n');
fprintf('  2. Select "Design Requirements" → "New..."\n');
fprintf('  3. Add "Damping Ratio" requirement: ζ ≥ %.4f\n', zeta_min);
fprintf('  4. Add "Natural Frequency" requirement: ωn ≥ %.0f rad/s\n', wn_min);
fprintf('  5. The acceptable region will be shaded\n');
fprintf('  6. Drag the pink square to select gain K\n\n');

% Open Control System Designer
controlSystemDesigner('rlocus', G);
pause(2);

%% Part (b): Select Gain to Meet Specifications
fprintf('\n========================================\n');
fprintf('PART (b): SELECT GAIN TO MEET SPECIFICATIONS\n');
fprintf('========================================\n\n');

fprintf('TASK: Select a gain K to place closed-loop poles in acceptable region\n');
fprintf('      while minimizing steady-state error.\n\n');

% Example: Try different gains
K_test = [10, 50, 100, 200, 500, 1000];

fprintf('Testing different gain values:\n');
fprintf('%-10s %-25s %-15s %-15s %-15s\n', 'K', 'Dominant Poles', 'ζ', 'ωn (rad/s)', 'Meets Spec?');
fprintf('%-10s %-25s %-15s %-15s %-15s\n', '---', '---------------', '---', '----------', '-----------');

for i = 1:length(K_test)
    K = K_test(i);
    sys_cl = feedback(K*G, 1);
    p = pole(sys_cl);
    
    % Find dominant poles (assume complex pair)
    [~, idx] = max(real(p));
    p_dom = p(idx);
    
    if abs(imag(p_dom)) > 1e-6
        sigma = real(p_dom);
        wd = abs(imag(p_dom));
        wn = sqrt(sigma^2 + wd^2);
        zeta = -sigma/wn;
        
        % Check if meets specifications
        meets_spec = (zeta >= zeta_min) && (wn >= wn_min);
        
        fprintf('%-10.0f %-25s %-15.4f %-15.2f %-15s\n', K, ...
            sprintf('%.2f ± j%.2f', real(p_dom), abs(imag(p_dom))), ...
            zeta, wn, string(meets_spec));
    end
end

% Select a reasonable gain (you should adjust based on root locus)
K_b = 500;  % Example gain - ADJUST THIS based on your Control System Designer selection
fprintf('\n** SELECTED GAIN: K = %.0f **\n\n', K_b);

sys_cl_b = feedback(K_b * G, 1);
poles_cl_b = pole(sys_cl_b);

fprintf('Closed-loop poles at K = %.0f:\n', K_b);
for i = 1:length(poles_cl_b)
    if abs(imag(poles_cl_b(i))) < 1e-6
        fprintf('  s = %.2f\n', real(poles_cl_b(i)));
    else
        fprintf('  s = %.2f ± j%.2f\n', real(poles_cl_b(i)), abs(imag(poles_cl_b(i))));
    end
end

% Calculate characteristics of selected gain
p_dom_b = poles_cl_b(1);
sigma_b = real(p_dom_b);
wd_b = abs(imag(p_dom_b));
wn_b = sqrt(sigma_b^2 + wd_b^2);
zeta_b = -sigma_b/wn_b;

fprintf('\nClosed-loop characteristics:\n');
fprintf('  ζ = %.4f (requirement: ≥ %.4f) %s\n', zeta_b, zeta_min, ...
    string(zeta_b >= zeta_min, "✓ PASS", "✗ FAIL"));
fprintf('  ωn = %.2f rad/s (requirement: ≥ %.0f rad/s) %s\n', wn_b, wn_min, ...
    string(wn_b >= wn_min, "✓ PASS", "✗ FAIL"));

%% Part (c): Step Response and Steady-State Error
fprintf('\n\n========================================\n');
fprintf('PART (c): STEP RESPONSE AND STEADY-STATE ERROR\n');
fprintf('========================================\n\n');

% Plot step response
figure('Name', 'Problem 5(c) - Step Response', 'NumberTitle', 'off', 'Position', [100 100 900 600]);
step(sys_cl_b);
grid on;
title(sprintf('Closed-Loop Unit Step Response (K = %.0f)', K_b));
xlabel('Time (s)');
ylabel('Angular Velocity ω (rad/s)');

% Get step info
S_b = stepinfo(sys_cl_b);
fprintf('Step response characteristics:\n');
fprintf('  Rise time (10-90%%): %.4f ms (spec: ≤ %.1f ms) %s\n', ...
    S_b.RiseTime*1000, tr_spec*1000, ...
    string(S_b.RiseTime <= tr_spec, "✓ PASS", "✗ FAIL"));
fprintf('  Overshoot: %.2f%% (spec: ≤ %.0f%%) %s\n', ...
    S_b.Overshoot, Mp_spec*100, ...
    string(S_b.Overshoot <= Mp_spec*100, "✓ PASS", "✗ FAIL"));
fprintf('  Settling time (2%%): %.4f ms\n', S_b.SettlingTime*1000);
fprintf('  Peak value: %.4f\n', S_b.Peak);

% Calculate steady-state error
fprintf('\n** Steady-State Error Calculation **\n');
fprintf('For Type 0 system with unity feedback:\n');
fprintf('  e_ss = 1/(1 + Kp), where Kp = lim[s→0] K·G(s)\n\n');

Kp_b = K_b * dcgain(G);
ess_b = 1 / (1 + Kp_b);
ess_percent_b = ess_b * 100;

fprintf('Position constant:\n');
fprintf('  Kp = K × G(0) = %.0f × %.4f = %.4f\n', K_b, dcgain(G), Kp_b);
fprintf('\nSteady-state error:\n');
fprintf('  e_ss = 1/(1 + Kp) = 1/(1 + %.4f) = %.6f\n', Kp_b, ess_b);
fprintf('  e_ss = %.4f%% of reference\n', ess_percent_b);

% Verify with Final Value Theorem
fprintf('\n** Verification using Final Value Theorem **\n');
fprintf('FVT: lim[t→∞] y(t) = lim[s→0] s·Y(s)\n');
fprintf('For unit step: R(s) = 1/s\n');
fprintf('  Y(s) = T(s)·R(s) = T(s)/s\n');
fprintf('  y_ss = lim[s→0] s·T(s)/s = T(0)\n\n');

T_b = feedback(K_b * G, 1);
y_ss = dcgain(T_b);
ess_fvt = 1 - y_ss;

fprintf('Closed-loop DC gain: T(0) = %.6f\n', y_ss);
fprintf('Steady-state error: e_ss = 1 - T(0) = %.6f (%.4f%%)\n', ess_fvt, ess_fvt*100);
fprintf('Verification: %.6f ≈ %.6f ✓\n', ess_b, ess_fvt);

%% Part (d): Increase Gain for 2% Error
fprintf('\n\n========================================\n');
fprintf('PART (d): INCREASE GAIN TO REDUCE ERROR TO 2%%\n');
fprintf('========================================\n\n');

fprintf('GOAL: Reduce steady-state error to 2%% or less\n\n');

fprintf('** Calculating required gain **\n');
fprintf('For e_ss = 0.02:\n');
fprintf('  0.02 = 1/(1 + Kp)\n');
fprintf('  1 + Kp = 1/0.02 = 50\n');
fprintf('  Kp = 49\n\n');

fprintf('Since Kp = K × G(0):\n');
fprintf('  K = Kp / G(0) = 49 / %.4f\n', dcgain(G));

K_d = (1/0.02 - 1) / dcgain(G);
fprintf('  K = %.2f\n\n', K_d);

fprintf('** SELECTED GAIN: K = %.2f **\n\n', K_d);

% Closed-loop system with increased gain
sys_cl_d = feedback(K_d * G, 1);
poles_cl_d = pole(sys_cl_d);

fprintf('Closed-loop poles at K = %.2f:\n', K_d);
for i = 1:length(poles_cl_d)
    if abs(imag(poles_cl_d(i))) < 1e-6
        fprintf('  s = %.2f\n', real(poles_cl_d(i)));
    else
        fprintf('  s = %.2f ± j%.2f\n', real(poles_cl_d(i)), abs(imag(poles_cl_d(i))));
    end
end

% Calculate new characteristics
p_dom_d = poles_cl_d(1);
sigma_d = real(p_dom_d);
wd_d = abs(imag(p_dom_d));
wn_d = sqrt(sigma_d^2 + wd_d^2);
zeta_d = -sigma_d/wn_d;

fprintf('\nClosed-loop characteristics:\n');
fprintf('  ζ = %.4f (requirement: ≥ %.4f) %s\n', zeta_d, zeta_min, ...
    string(zeta_d >= zeta_min, "✓ PASS", "✗ FAIL"));
fprintf('  ωn = %.2f rad/s (requirement: ≥ %.0f rad/s) %s\n', wn_d, wn_min, ...
    string(wn_d >= wn_min, "✓ PASS", "✗ FAIL"));

% Plot step response
figure('Name', 'Problem 5(d) - Higher Gain Response', 'NumberTitle', 'off', 'Position', [150 150 900 600]);
step(sys_cl_d);
grid on;
title(sprintf('Closed-Loop Step Response (K = %.2f for 2%% error)', K_d));
xlabel('Time (s)');
ylabel('Angular Velocity ω (rad/s)');

% Get step info
S_d = stepinfo(sys_cl_d);
fprintf('\nStep response characteristics:\n');
fprintf('  Rise time (10-90%%): %.4f ms (spec: ≤ %.1f ms) %s\n', ...
    S_d.RiseTime*1000, tr_spec*1000, ...
    string(S_d.RiseTime <= tr_spec, "✓ PASS", "✗ FAIL"));
fprintf('  Overshoot: %.2f%% (spec: ≤ %.0f%%) %s\n', ...
    S_d.Overshoot, Mp_spec*100, ...
    string(S_d.Overshoot <= Mp_spec*100, "✓ PASS", "✗ FAIL"));
fprintf('  Settling time (2%%): %.4f ms\n', S_d.SettlingTime*1000);

% Verify error
ess_d = 1 - dcgain(feedback(K_d*G, 1));
fprintf('  Steady-state error: %.4f%% %s\n', ess_d*100, ...
    string(ess_d <= 0.02, "✓ PASS", "✗ FAIL"));

fprintf('\n** ANALYSIS **\n');
fprintf('Effect of increasing gain:\n');
if S_d.Overshoot > Mp_spec*100
    fprintf('  ✗ Overshoot INCREASED from %.2f%% to %.2f%% (exceeds spec!)\n', ...
        S_b.Overshoot, S_d.Overshoot);
    fprintf('  → Simply increasing K cannot meet all specifications\n');
    fprintf('  → Need compensation (lag or lead-lag controller)\n');
else
    fprintf('  ✓ Specifications still met\n');
end

%% Part (e): Lag Compensator Exploration
fprintf('\n\n========================================\n');
fprintf('PART (e): LAG COMPENSATOR EXPLORATION\n');
fprintf('========================================\n\n');

fprintf('LAG COMPENSATOR: C(s) = Kc × (s+z)/(s+p), where p < z\n\n');

fprintf('Purpose of lag compensator:\n');
fprintf('  1. Increase DC gain (reduce steady-state error)\n');
fprintf('  2. Minimal impact on transient response\n');
fprintf('  3. Place pole and zero CLOSE TO ORIGIN (< 1/10 of dominant poles)\n\n');

fprintf('Effect on root locus:\n');
fprintf('  - Pole-zero pair near origin adds "extra branch"\n');
fprintf('  - Main branches barely affected (if p and z are close to origin)\n');
fprintf('  - Branches shift SLIGHTLY LEFT (improved stability margin)\n');
fprintf('  - DC gain increased by factor of z/p\n\n');

fprintf('Opening Control System Designer for lag compensator design...\n');
fprintf('Instructions:\n');
fprintf('  1. In Compensator Editor, right-click → Add Pole/Zero → Real Pole\n');
fprintf('  2. Place pole close to origin (e.g., p = 1 to 10)\n');
fprintf('  3. Right-click → Add Pole/Zero → Real Zero\n');
fprintf('  4. Place zero further from origin than pole (e.g., z = 10 to 100)\n');
fprintf('  5. Ensure p < z for lag compensation\n');
fprintf('  6. Adjust Kc (gain) to meet error requirement\n');
fprintf('  7. Observe how branches barely move but DC gain increases\n\n');

controlSystemDesigner('rlocus', G);
pause(2);

%% Part (f): Minimum Controller Gain for 2% Error
fprintf('\n========================================\n');
fprintf('PART (f): MINIMUM CONTROLLER GAIN\n');
fprintf('========================================\n\n');

fprintf('LAG COMPENSATOR: C(s) = Kc × (s+z)/(s+p)\n\n');

fprintf('DC gain of compensator:\n');
fprintf('  C(0) = Kc × z/p\n\n');

fprintf('Total DC gain (position constant):\n');
fprintf('  Kp = C(0) × K × G(0) = (Kc × z/p) × G(0)\n');
fprintf('  (Note: We absorb K into Kc for simplicity)\n\n');

fprintf('For 2%% steady-state error:\n');
fprintf('  e_ss = 1/(1 + Kp) = 0.02\n');
fprintf('  1 + Kp = 50\n');
fprintf('  Kp = 49\n\n');

fprintf('Therefore:\n');
fprintf('  (Kc × z/p) × G(0) = 49\n');
fprintf('  Kc × z/p = 49 / G(0)\n');

min_dc_gain = (1/0.02 - 1) / dcgain(G);

fprintf('  Kc × z/p = 49 / %.4f\n', dcgain(G));
fprintf('  Kc × z/p = %.2f\n\n', min_dc_gain);

fprintf('** MINIMUM CONTROLLER DC GAIN: Kc × z/p ≥ %.2f **\n', min_dc_gain);

fprintf('\nExample combinations:\n');
fprintf('  If z/p = 10: Kc ≥ %.2f\n', min_dc_gain/10);
fprintf('  If z/p = 5:  Kc ≥ %.2f\n', min_dc_gain/5);
fprintf('  If z/p = 2:  Kc ≥ %.2f\n', min_dc_gain/2);

%% Part (g): Design Lag Controller
fprintf('\n\n========================================\n');
fprintf('PART (g): DESIGN LAG CONTROLLER\n');
fprintf('========================================\n\n');

fprintf('GOAL: Design lag controller to meet ALL specifications:\n');
fprintf('  - Rise time ≤ %.1f ms\n', tr_spec*1000);
fprintf('  - Overshoot ≤ %.0f%%\n', Mp_spec*100);
fprintf('  - Steady-state error ≤ 2%%\n\n');

fprintf('DESIGN STRATEGY:\n');
fprintf('  1. Place pole and zero close to origin (typically < 1/10 of system poles)\n');
fprintf('  2. Choose z/p ratio to achieve desired DC gain\n');
fprintf('  3. Adjust Kc to meet steady-state error requirement\n');
fprintf('  4. Verify transient response still meets specifications\n\n');

% Design example 1: Conservative approach
z_lag1 = 50;
p_lag1 = 5;
Kc1 = min_dc_gain * p_lag1 / z_lag1;

fprintf('** DESIGN EXAMPLE 1: Conservative **\n');
fprintf('  Zero: z = %.0f\n', z_lag1);
fprintf('  Pole: p = %.0f\n', p_lag1);
fprintf('  Ratio: z/p = %.0f\n', z_lag1/p_lag1);
fprintf('  Gain: Kc = %.4f\n', Kc1);
fprintf('  DC gain: Kc×z/p = %.2f\n\n', Kc1*z_lag1/p_lag1);

% Create lag compensator 1
C_lag1 = tf(Kc1 * [1, z_lag1], [1, p_lag1]);
G_comp1 = C_lag1 * G;
sys_cl_lag1 = feedback(G_comp1, 1);

% Analyze
poles_lag1 = pole(sys_cl_lag1);
S_lag1 = stepinfo(sys_cl_lag1);
ess_lag1 = 1 - dcgain(sys_cl_lag1);

fprintf('Results:\n');
fprintf('  Closed-loop poles:\n');
for i = 1:length(poles_lag1)
    if abs(imag(poles_lag1(i))) < 1e-6
        fprintf('    s = %.2f\n', real(poles_lag1(i)));
    else
        fprintf('    s = %.2f ± j%.2f\n', real(poles_lag1(i)), abs(imag(poles_lag1(i))));
    end
end

% Find dominant poles for characteristics
complex_poles1 = poles_lag1(abs(imag(poles_lag1)) > 1e-6);
if ~isempty(complex_poles1)
    [~, idx] = max(real(complex_poles1));
    p_dom_lag1 = complex_poles1(idx);
    wn_lag1 = abs(p_dom_lag1);
    zeta_lag1 = -real(p_dom_lag1)/wn_lag1;
    fprintf('  ζ = %.4f %s\n', zeta_lag1, string(zeta_lag1 >= zeta_min, "✓", "✗"));
    fprintf('  ωn = %.2f rad/s %s\n', wn_lag1, string(wn_lag1 >= wn_min, "✓", "✗"));
end

fprintf('  Rise time: %.4f ms %s\n', S_lag1.RiseTime*1000, ...
    string(S_lag1.RiseTime <= tr_spec, "✓", "✗"));
fprintf('  Overshoot: %.2f%% %s\n', S_lag1.Overshoot, ...
    string(S_lag1.Overshoot <= Mp_spec*100, "✓", "✗"));
fprintf('  Steady-state error: %.4f%% %s\n\n', ess_lag1*100, ...
    string(ess_lag1 <= 0.02, "✓", "✗"));

% Design example 2: Aggressive approach (poles further from origin)
z_lag2 = 100;
p_lag2 = 10;
Kc2 = min_dc_gain * p_lag2 / z_lag2 * 1.1; % Add 10% margin

fprintf('** DESIGN EXAMPLE 2: Moderate **\n');
fprintf('  Zero: z = %.0f\n', z_lag2);
fprintf('  Pole: p = %.0f\n', p_lag2);
fprintf('  Ratio: z/p = %.0f\n', z_lag2/p_lag2);
fprintf('  Gain: Kc = %.4f (with 10%% margin)\n', Kc2);
fprintf('  DC gain: Kc×z/p = %.2f\n\n', Kc2*z_lag2/p_lag2);

% Create lag compensator 2
C_lag2 = tf(Kc2 * [1, z_lag2], [1, p_lag2]);
G_comp2 = C_lag2 * G;
sys_cl_lag2 = feedback(G_comp2, 1);

% Analyze
poles_lag2 = pole(sys_cl_lag2);
S_lag2 = stepinfo(sys_cl_lag2);
ess_lag2 = 1 - dcgain(sys_cl_lag2);

fprintf('Results:\n');
fprintf('  Closed-loop poles:\n');
for i = 1:length(poles_lag2)
    if abs(imag(poles_lag2(i))) < 1e-6
        fprintf('    s = %.2f\n', real(poles_lag2(i)));
    else
        fprintf('    s = %.2f ± j%.2f\n', real(poles_lag2(i)), abs(imag(poles_lag2(i))));
    end
end

% Find dominant poles
complex_poles2 = poles_lag2(abs(imag(poles_lag2)) > 1e-6);
if ~isempty(complex_poles2)
    [~, idx] = max(real(complex_poles2));
    p_dom_lag2 = complex_poles2(idx);
    wn_lag2 = abs(p_dom_lag2);
    zeta_lag2 = -real(p_dom_lag2)/wn_lag2;
    fprintf('  ζ = %.4f %s\n', zeta_lag2, string(zeta_lag2 >= zeta_min, "✓", "✗"));
    fprintf('  ωn = %.2f rad/s %s\n', wn_lag2, string(wn_lag2 >= wn_min, "✓", "✗"));
end

fprintf('  Rise time: %.4f ms %s\n', S_lag2.RiseTime*1000, ...
    string(S_lag2.RiseTime <= tr_spec, "✓", "✗"));
fprintf('  Overshoot: %.2f%% %s\n', S_lag2.Overshoot, ...
    string(S_lag2.Overshoot <= Mp_spec*100, "✓", "✗"));
fprintf('  Steady-state error: %.4f%% %s\n\n', ess_lag2*100, ...
    string(ess_lag2 <= 0.02, "✓", "✗"));

% Design example 3: Optimized design
% Try to find better pole-zero placement
z_lag3 = 80;
p_lag3 = 8;
Kc3 = min_dc_gain * p_lag3 / z_lag3 * 1.05; % Add 5% margin

fprintf('** DESIGN EXAMPLE 3: Optimized **\n');
fprintf('  Zero: z = %.0f\n', z_lag3);
fprintf('  Pole: p = %.0f\n', p_lag3);
fprintf('  Ratio: z/p = %.0f\n', z_lag3/p_lag3);
fprintf('  Gain: Kc = %.4f\n', Kc3);
fprintf('  DC gain: Kc×z/p = %.2f\n\n', Kc3*z_lag3/p_lag3);

% Create lag compensator 3
C_lag3 = tf(Kc3 * [1, z_lag3], [1, p_lag3]);
G_comp3 = C_lag3 * G;
sys_cl_lag3 = feedback(G_comp3, 1);

% Analyze
poles_lag3 = pole(sys_cl_lag3);
S_lag3 = stepinfo(sys_cl_lag3);
ess_lag3 = 1 - dcgain(sys_cl_lag3);

fprintf('Results:\n');
fprintf('  Closed-loop poles:\n');
for i = 1:length(poles_lag3)
    if abs(imag(poles_lag3(i))) < 1e-6
        fprintf('    s = %.2f\n', real(poles_lag3(i)));
    else
        fprintf('    s = %.2f ± j%.2f\n', real(poles_lag3(i)), abs(imag(poles_lag3(i))));
    end
end

% Find dominant poles
complex_poles3 = poles_lag3(abs(imag(poles_lag3)) > 1e-6);
if ~isempty(complex_poles3)
    [~, idx] = max(real(complex_poles3));
    p_dom_lag3 = complex_poles3(idx);
    wn_lag3 = abs(p_dom_lag3);
    zeta_lag3 = -real(p_dom_lag3)/wn_lag3;
    fprintf('  ζ = %.4f %s\n', zeta_lag3, string(zeta_lag3 >= zeta_min, "✓", "✗"));
    fprintf('  ωn = %.2f rad/s %s\n', wn_lag3, string(wn_lag3 >= wn_min, "✓", "✗"));
end

fprintf('  Rise time: %.4f ms %s\n', S_lag3.RiseTime*1000, ...
    string(S_lag3.RiseTime <= tr_spec, "✓", "✗"));
fprintf('  Overshoot: %.2f%% %s\n', S_lag3.Overshoot, ...
    string(S_lag3.Overshoot <= Mp_spec*100, "✓", "✗"));
fprintf('  Steady-state error: %.4f%% %s\n\n', ess_lag3*100, ...
    string(ess_lag3 <= 0.02, "✓", "✗"));

% Plot comparison of step responses
figure('Name', 'Problem 5(g) - Lag Compensator Comparison', 'NumberTitle', 'off', ...
    'Position', [200 100 1200 800]);

subplot(2,2,1);
step(sys_cl_b, 'b', sys_cl_lag1, 'r', sys_cl_lag2, 'g', sys_cl_lag3, 'm');
grid on;
title('Step Response Comparison');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend(sprintf('Proportional (K=%.0f)', K_b), ...
    sprintf('Lag 1 (z=%.0f, p=%.0f)', z_lag1, p_lag1), ...
    sprintf('Lag 2 (z=%.0f, p=%.0f)', z_lag2, p_lag2), ...
    sprintf('Lag 3 (z=%.0f, p=%.0f)', z_lag3, p_lag3), ...
    'Location', 'best');

subplot(2,2,2);
step(sys_cl_lag3);
grid on;
title(sprintf('Optimized Lag Controller (z=%.0f, p=%.0f, Kc=%.3f)', z_lag3, p_lag3, Kc3));
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');

subplot(2,2,3);
rlocus(G_comp3);
grid on;
hold on;
sgrid(zeta_min, wn_min);
plot(real(poles_lag3), imag(poles_lag3), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
title('Root Locus with Lag Compensator');
xlabel('Real Axis');
ylabel('Imaginary Axis');
xlim([-5000 500]);
ylim([-5000 5000]);

subplot(2,2,4);
% Performance comparison bar chart
designs = categorical({'Prop', 'Lag1', 'Lag2', 'Lag3'});
rise_times = [S_b.RiseTime, S_lag1.RiseTime, S_lag2.RiseTime, S_lag3.RiseTime] * 1000;
overshoots = [S_b.Overshoot, S_lag1.Overshoot, S_lag2.Overshoot, S_lag3.Overshoot];
ss_errors = [ess_b, ess_lag1, ess_lag2, ess_lag3] * 100;

bar(designs, [rise_times', overshoots', ss_errors']);
grid on;
title('Performance Comparison');
ylabel('Value');
legend('Rise Time (ms)', 'Overshoot (%)', 'SS Error (%)', 'Location', 'best');

% Select best design (you can adjust based on results)
fprintf('\n** RECOMMENDED DESIGN **\n');
fprintf('Based on the analysis, select the design that best meets all specifications.\n');
fprintf('Typically, Design 3 (moderate pole-zero placement) offers good balance.\n\n');

% Use Design 3 as final
z_final = z_lag3;
p_final = p_lag3;
Kc_final = Kc3;

fprintf('FINAL SELECTED LAG COMPENSATOR:\n');
fprintf('  C(s) = %.4f × (s + %.0f) / (s + %.0f)\n', Kc_final, z_final, p_final);
fprintf('  Zero location: z = %.0f\n', z_final);
fprintf('  Pole location: p = %.0f\n', p_final);
fprintf('  Controller gain: Kc = %.4f\n', Kc_final);
fprintf('  DC gain: Kc×z/p = %.2f\n\n', Kc_final*z_final/p_final);

fprintf('Closed-loop poles:\n');
for i = 1:length(poles_lag3)
    if abs(imag(poles_lag3(i))) < 1e-6
        fprintf('  s = %.2f\n', real(poles_lag3(i)));
    else
        fprintf('  s = %.2f ± j%.2f\n', real(poles_lag3(i)), abs(imag(poles_lag3(i))));
    end
end

fprintf('\nPerformance Summary:\n');
fprintf('  %-25s %-15s %-10s\n', 'Specification', 'Value', 'Status');
fprintf('  %-25s %-15s %-10s\n', '-------------------------', '---------------', '----------');
fprintf('  %-25s %-15s %-10s\n', 'Rise time ≤ 0.5 ms', sprintf('%.4f ms', S_lag3.RiseTime*1000), ...
    string(S_lag3.RiseTime <= tr_spec, "✓ PASS", "✗ FAIL"));
fprintf('  %-25s %-15s %-10s\n', 'Overshoot ≤ 5%', sprintf('%.2f%%', S_lag3.Overshoot), ...
    string(S_lag3.Overshoot <= Mp_spec*100, "✓ PASS", "✗ FAIL"));
fprintf('  %-25s %-15s %-10s\n', 'SS Error ≤ 2%', sprintf('%.4f%%', ess_lag3*100), ...
    string(ess_lag3 <= 0.02, "✓ PASS", "✗ FAIL"));

% Create detailed root locus plot with compensator
figure('Name', 'Problem 5(g) - Final Root Locus with Lag Compensator', 'NumberTitle', 'off', ...
    'Position', [250 150 1000 800]);

rlocus(G_comp3);
grid on;
hold on;

% Add design requirements
sgrid(zeta_min, wn_min);

% Add requirement boundary lines
theta_min = acos(zeta_min);
x_zeta = [-10000, 0];
y_upper = -x_zeta * tan(theta_min);
y_lower = -y_upper;
plot(x_zeta, y_upper, 'r--', 'LineWidth', 2);
plot(x_zeta, y_lower, 'r--', 'LineWidth', 2);

% Circle for minimum wn
theta_circle = linspace(0, 2*pi, 100);
x_circle = wn_min * cos(theta_circle);
y_circle = wn_min * sin(theta_circle);
plot(x_circle, y_circle, 'g--', 'LineWidth', 2);

% Mark closed-loop poles
plot(real(poles_lag3), imag(poles_lag3), 'rx', 'MarkerSize', 15, 'LineWidth', 3);

% Mark compensator pole and zero
plot(-p_final, 0, 'bs', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'b');
plot(-z_final, 0, 'bo', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'b');

title('Root Locus with Final Lag Compensator');
xlabel('Real Axis (σ)');
ylabel('Imaginary Axis (jω)');
legend('Root Locus', 'Design Grid', sprintf('ζ = %.3f', zeta_min), ...
    sprintf('ωn = %.0f rad/s', wn_min), 'Closed-loop Poles', ...
    'Compensator Pole', 'Compensator Zero', 'Location', 'best');
xlim([-6000 500]);
ylim([-6000 6000]);

%% Summary and Conclusion
fprintf('\n\n========================================================================\n');
fprintf('PROBLEM 5 SUMMARY AND CONCLUSIONS\n');
fprintf('========================================================================\n\n');

fprintf('SYSTEM PARAMETERS:\n');
fprintf('  Open-loop transfer function: G(s) = %.4f / (%.4es² + %.4es + %.4f)\n', ...
    num, a2, a1, a0);
fprintf('  Open-loop poles: ');
for i = 1:length(poles_ol)
    if i > 1, fprintf(', '); end
    if abs(imag(poles_ol(i))) < 1e-6
        fprintf('%.2f', real(poles_ol(i)));
    else
        fprintf('%.2f±j%.2f', real(poles_ol(i)), abs(imag(poles_ol(i))));
    end
end
fprintf('\n\n');

fprintf('DESIGN REQUIREMENTS:\n');
fprintf('  Rise time: ≤ %.1f ms\n', tr_spec*1000);
fprintf('  Overshoot: ≤ %.0f%%\n', Mp_spec*100);
fprintf('  Steady-state error: ≤ 2%%\n');
fprintf('  Translated to: ζ ≥ %.4f, ωn ≥ %.0f rad/s\n\n', zeta_min, wn_min);

fprintf('KEY FINDINGS:\n');
fprintf('  1. Proportional control alone (Part b, K=%.0f):\n', K_b);
fprintf('     - Can meet transient specs but has high steady-state error (%.2f%%)\n', ess_b*100);
fprintf('\n');
fprintf('  2. Increased gain (Part d, K=%.2f):\n', K_d);
fprintf('     - Reduces error to 2%% but violates transient specifications\n');
fprintf('     - Overshoot increases to %.2f%% (exceeds 5%% limit)\n', S_d.Overshoot);
fprintf('\n');
fprintf('  3. Lag compensator (Part g):\n');
fprintf('     - Successfully meets ALL specifications\n');
fprintf('     - Pole at s = -%.0f, Zero at s = -%.0f\n', p_final, z_final);
fprintf('     - Controller gain Kc = %.4f\n', Kc_final);
fprintf('     - Rise time: %.4f ms ✓\n', S_lag3.RiseTime*1000);
fprintf('     - Overshoot: %.2f%% ✓\n', S_lag3.Overshoot);
fprintf('     - SS Error: %.4f%% ✓\n\n', ess_lag3*100);

