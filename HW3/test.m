%% ME351 HW3 - Complete Solutions for Problems 1-5
% Root Locus Analysis and Design

clear all; close all; clc;

%% ========================================================================
%% PROBLEM 2
%% ========================================================================
fprintf('\n========================================\n');
fprintf('PROBLEM 2\n');
fprintf('========================================\n');

% TODO: Enter your transfer function from the problem image
% Example: G2(s) = K / ((s+1)(s+3)(s+5))
num2 = 1;
den2 = conv([1 1], conv([1 3], [1 5]));
G2 = tf(num2, den2);

fprintf('Transfer Function:\n');
disp(G2);

%% Part (a): Root Locus Sketch
fprintf('\n--- Part (a): Root Locus Sketch ---\n');

poles2 = pole(G2);
zeros2 = zero(G2);
n_poles2 = length(poles2);
n_zeros2 = length(zeros2);

fprintf('Poles: ');
disp(poles2');
fprintf('Zeros: ');
if isempty(zeros2)
    fprintf('None\n');
else
    disp(zeros2');
end

% Real-axis segments
fprintf('\nReal-axis segments: Root locus to the left of ODD number of singularities\n');

% Asymptotes
n_asymptotes2 = n_poles2 - n_zeros2;
asymptote_angles2 = zeros(1, n_asymptotes2);
for k = 1:n_asymptotes2
    asymptote_angles2(k) = (2*k - 1) * 180 / n_asymptotes2;
end
fprintf('Asymptote angles: ');
fprintf('%.1f° ', asymptote_angles2);
fprintf('\n');

centroid2 = (sum(poles2) - sum(zeros2)) / (n_poles2 - n_zeros2);
fprintf('Centroid: σ = %.4f\n', real(centroid2));

% Breakaway points
fprintf('\nBreakaway/break-in points: Solve dK/ds = 0\n');

%% Part (b): Stability Range using Routh Array
fprintf('\n--- Part (b): Stability Range (Routh Array) ---\n');

% Characteristic equation
char_eq2 = den2;
fprintf('Characteristic polynomial: ');
for i = 1:length(char_eq2)
    if i < length(char_eq2)
        fprintf('%.4f*s^%d + ', char_eq2(i), length(char_eq2)-i);
    else
        fprintf('%.4f + K*%.4f = 0\n', char_eq2(i), num2(end));
    end
end

% Build Routh array
fprintf('\nBuilding Routh array...\n');
n = length(char_eq2);
routh = zeros(n, ceil(n/2));

% First two rows
routh(1, :) = [char_eq2(1:2:end), zeros(1, ceil(n/2)-length(char_eq2(1:2:end)))];
routh(2, :) = [char_eq2(2:2:end), zeros(1, ceil(n/2)-length(char_eq2(2:2:end)))];

fprintf('Routh Array (without K in last term):\n');
fprintf('s^%d: ', n-1);
fprintf('%.4f  ', routh(1,:));
fprintf('\n');
fprintf('s^%d: ', n-2);
fprintf('%.4f  ', routh(2,:));
fprintf('\n');

fprintf('\nFor stability, all elements in first column must be POSITIVE\n');
fprintf('Find range of K that makes all first column elements > 0\n');
fprintf('System is UNSTABLE when K is outside this range\n');

%% Part (c): Gain for fastest settling time
fprintf('\n--- Part (c): Fastest Settling Time ---\n');
fprintf('Settling time ts ≈ 4/|Re(s)|, where s is the dominant pole\n');
fprintf('For fastest settling time, move dominant poles as far LEFT as possible\n');
fprintf('Find the K value that maximizes |Re(s)| of the dominant pole\n');
fprintf('This typically occurs just before instability or at a local maximum\n');

% MATLAB Root Locus Plot
figure('Name', 'Problem 2 - Root Locus', 'NumberTitle', 'off');
rlocus(G2);
grid on;
title('Problem 2: Root Locus');
xlabel('Real Axis');
ylabel('Imaginary Axis');
sgrid;

% Open Control System Designer
controlSystemDesigner('rlocus', G2);
fprintf('\nControl System Designer opened for Problem 2\n');

%% ========================================================================
%% PROBLEM 3
%% ========================================================================
fprintf('\n========================================\n');
fprintf('PROBLEM 3\n');
fprintf('========================================\n');

% TODO: Enter your transfer function from the problem image
% This should have (s+2) in numerator compared to previous problems
% Example: G3(s) = K*(s+2) / (s(s+4)(s+6))
num3 = [1 2];  % (s+2)
den3 = conv([1 0], conv([1 4], [1 6]));
G3 = tf(num3, den3);

fprintf('Transfer Function:\n');
disp(G3);

%% Part (a): Root Locus Sketch
fprintf('\n--- Part (a): Root Locus Sketch ---\n');

poles3 = pole(G3);
zeros3 = zero(G3);

fprintf('Poles: ');
disp(poles3');
fprintf('Zeros: ');
disp(zeros3');

% Complete root locus analysis
n_poles3 = length(poles3);
n_zeros3 = length(zeros3);

% Asymptotes
n_asymptotes3 = n_poles3 - n_zeros3;
asymptote_angles3 = zeros(1, n_asymptotes3);
for k = 1:n_asymptotes3
    asymptote_angles3(k) = (2*k - 1) * 180 / n_asymptotes3;
end
fprintf('\nAsymptote angles: ');
fprintf('%.1f° ', asymptote_angles3);
fprintf('\n');

centroid3 = (sum(poles3) - sum(zeros3)) / (n_poles3 - n_zeros3);
fprintf('Centroid: σ = %.4f\n', real(centroid3));

% MATLAB Root Locus Plot
figure('Name', 'Problem 3 - Root Locus', 'NumberTitle', 'off');
rlocus(G3);
grid on;
title('Problem 3: Root Locus (with zero at s=-2)');
xlabel('Real Axis');
ylabel('Imaginary Axis');
sgrid;

% Open Control System Designer
controlSystemDesigner('rlocus', G3);
fprintf('\nControl System Designer opened for Problem 3\n');

%% Part (b): Compare root loci from Problems 1-3
fprintf('\n--- Part (b): Comparison of Root Loci (Effect of (s+2) zero) ---\n');
fprintf('\nComparing Problems 1-3:\n');
fprintf('- The zero in the numerator (s+2) pulls the root locus branches\n');
fprintf('- Zeros attract root locus branches toward them\n');
fprintf('- A zero on the real axis typically:\n');
fprintf('  1. Reduces number of asymptotes\n');
fprintf('  2. Pulls branches to the LEFT (improves stability/settling time)\n');
fprintf('  3. Can eliminate some breakaway points\n');
fprintf('  4. May prevent branches from going into RHP\n');
fprintf('\nKey observation: Adding a LHP zero generally IMPROVES system stability\n');

% Overlay comparison plot
figure('Name', 'Comparison - Problems 1-3', 'NumberTitle', 'off');
rlocus(G1, 'b', G2, 'r', G3, 'g');
grid on;
title('Comparison: Effect of Zero (s+2)');
xlabel('Real Axis');
ylabel('Imaginary Axis');
legend('Problem 1', 'Problem 2', 'Problem 3 (with zero)', 'Location', 'best');
sgrid;

%% ========================================================================
%% PROBLEM 4: Satellite Attitude Control
%% ========================================================================
fprintf('\n========================================\n');
fprintf('PROBLEM 4: Satellite Attitude Control\n');
fprintf('========================================\n');

% Given
J = 1;  % kg-m^2
G_sat = tf(1, [J, 0, 0]);  % G(s) = 1/(J*s^2)

fprintf('Plant transfer function: G(s) = 1/(s^2)\n');
disp(G_sat);

%% Part (a): Proportional Control (C(s) = Kp)
fprintf('\n--- Part (a): Proportional Control Only ---\n');

% For C(s) = Kp, open-loop is Kp*G(s) = Kp/s^2
G_prop = G_sat;

fprintf('Open-loop with P-control: Kp/s^2\n');
fprintf('Poles: s = 0 (double pole at origin)\n');
fprintf('Zeros: None\n');

% Root locus analysis
fprintf('\nRoot locus characteristics:\n');
fprintf('- 2 poles at origin, 0 zeros\n');
fprintf('- 2 asymptotes at ±90°\n');
fprintf('- Centroid: σ = 0/2 = 0\n');
fprintf('- Branches leave origin at ±90° (imaginary axis)\n');
fprintf('- Root locus stays on imaginary axis for all K\n');

figure('Name', 'Problem 4a - Proportional Control', 'NumberTitle', 'off');
rlocus(G_prop);
grid on;
title('Problem 4(a): Root Locus with Proportional Control (C(s) = K_p)');
xlabel('Real Axis');
ylabel('Imaginary Axis');
axis([-2 2 -3 3]);
sgrid;

% Open Control System Designer
controlSystemDesigner('rlocus', G_prop);
fprintf('\nControl System Designer opened for Problem 4(a)\n');

%% Part (b): Comment on response with varying Kp
fprintf('\n--- Part (b): Nature of Response ---\n');
fprintf('With proportional control only:\n');
fprintf('- All closed-loop poles remain on the imaginary axis (pure oscillation)\n');
fprintf('- System is MARGINALLY STABLE for all Kp > 0\n');
fprintf('- No damping - system oscillates indefinitely\n');
fprintf('- As Kp increases, oscillation frequency increases: ω = sqrt(Kp)\n');
fprintf('- Cannot achieve stable, damped response with P-control alone\n');

%% Part (c): PD Control characteristic equation
fprintf('\n--- Part (c): PD Control with Kp/Kd = 1 ---\n');

% C(s) = Kp + Kd*s, with Kp/Kd = 1, so Kp = Kd
% Let K = Kd, then C(s) = K + K*s = K(1 + s)
% Open-loop: K(1+s)/s^2 = K(s+1)/s^2

fprintf('C(s) = Kp + Kd*s, with K = Kd and Kp/Kd = 1\n');
fprintf('Therefore: C(s) = K(s + 1)\n');
fprintf('Open-loop: L(s) = C(s)*G(s) = K(s+1)/s^2\n');
fprintf('Characteristic equation: 1 + K(s+1)/s^2 = 0\n');
fprintf('Root locus form: s^2 + K(s+1) = 0\n');

G_pd = tf([1 1], [1 0 0]);  % (s+1)/s^2

%% Part (d): Root locus for PD control
fprintf('\n--- Part (d): Root Locus with PD Control ---\n');

poles_pd = [0; 0];
zeros_pd = -1;

fprintf('Poles: 0, 0 (double pole at origin)\n');
fprintf('Zeros: -1\n');
fprintf('Number of asymptotes: 2 - 1 = 1\n');
fprintf('Asymptote angle: 180°\n');
fprintf('Centroid: (0+0-(-1))/(2-1) = 1\n');

% Breakaway point
fprintf('\nBreakaway point: Solve dK/ds = 0\n');
fprintf('K = -s^2/(s+1)\n');
fprintf('dK/ds = -(2s(s+1) - s^2)/(s+1)^2 = -(s^2 + 2s)/(s+1)^2 = 0\n');
fprintf('s^2 + 2s = 0 → s(s+2) = 0\n');
fprintf('Breakaway points: s = 0 (at pole), s = -2\n');
K_breakaway_pd = -(-2)^2 / (-2+1);
fprintf('At s = -2: K = %.4f\n', K_breakaway_pd);

% Departure angle from double pole
fprintf('\nDeparture angles from double pole at origin:\n');
fprintf('θ_d = (180° + Σ∠zeros - Σ∠other poles) / (multiplicity)\n');
fprintf('θ_d = (180° + 180° - 0°) / 1 = 360° (or 0°) NOT CORRECT FOR DOUBLE POLE\n');
fprintf('For double pole: ±90° departure (perpendicular to real axis)\n');

figure('Name', 'Problem 4d - PD Control', 'NumberTitle', 'off');
rlocus(G_pd);
grid on;
title('Problem 4(d): Root Locus with PD Control (K = K_d, K_p/K_d = 1)');
xlabel('Real Axis');
ylabel('Imaginary Axis');
sgrid;

% Add annotations
hold on;
plot(-2, 0, 'ms', 'MarkerSize', 12, 'LineWidth', 2);
text(-2, 0.3, 'Breakaway point', 'FontSize', 10);

% Open Control System Designer
controlSystemDesigner('rlocus', G_pd);
fprintf('\nControl System Designer opened for Problem 4(d)\n');

%% Part (e): Effect of PD Control
fprintf('\n--- Part (e): Effect of PD Control ---\n');
fprintf('Comparing P-control vs PD-control:\n');
fprintf('- PD control adds a ZERO at s = -1\n');
fprintf('- This zero pulls the root locus branches to the LEFT\n');
fprintf('- Branches now move into the LEFT HALF PLANE\n');
fprintf('- System becomes STABLE with damping\n');
fprintf('- Can now achieve desired transient response (overshoot, settling time)\n');
fprintf('- The derivative term provides damping to the system\n');

%% Part (f): Lead Compensation
fprintf('\n--- Part (f): Lead Compensation ---\n');

% C(s) = Kc * (s+z)/(s+p) with z < p (for lead)
fprintf('Lead compensator: C(s) = Kc * (s+z)/(s+p)\n');
fprintf('Open-loop: L(s) = Kc * (s+z)/(s+p) * 1/s^2 = Kc(s+z)/(s^2(s+p))\n');
fprintf('Characteristic equation: 1 + Kc(s+z)/(s^2(s+p)) = 0\n');
fprintf('Root locus form: s^2(s+p) + Kc(s+z) = 0\n');
fprintf('Or: 1 + Kc(s+z)/(s^2(s+p)) = 0\n');

%% Part (g): Root loci for different pole locations
fprintf('\n--- Part (g): Lead Compensation with z=1, varying p ---\n');

z_lead = 1;
p_values = [20, 9, 3];

figure('Name', 'Problem 4g - Lead Compensation Comparison', 'NumberTitle', 'off');

for idx = 1:length(p_values)
    p_lead = p_values(idx);
    
    fprintf('\n(g.%d) p = %d:\n', idx, p_lead);
    
    % Transfer function
    G_lead = tf([1 z_lead], conv([1 0 0], [1 p_lead]));
    
    % Root locus analysis
    poles_lead = [0; 0; -p_lead];
    zeros_lead = -z_lead;
    
    fprintf('  Poles: 0, 0, %.0f\n', -p_lead);
    fprintf('  Zeros: %.0f\n', -z_lead);
    
    % Asymptotes
    n_asymp = 3 - 1;
    asymp_angles = [(2*1-1)*180/n_asymp, (2*2-1)*180/n_asymp];
    fprintf('  Asymptotes: %d branches at %.0f° and %.0f°\n', n_asymp, asymp_angles(1), asymp_angles(2));
    
    % Centroid
    cent = (sum(poles_lead) - sum(zeros_lead)) / n_asymp;
    fprintf('  Centroid: σ = (0+0+(%.0f)-(%.0f))/2 = %.2f\n', -p_lead, -z_lead, cent);
    
    % Breakaway points
    fprintf('  Breakaway points: Solve dK/ds = 0\n');
    fprintf('  (requires solving cubic equation)\n');
    
    % Plot
    subplot(1, 3, idx);
    rlocus(G_lead);
    grid on;
    title(sprintf('p = %d', p_lead));
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    sgrid;
    axis equal;
    xlim([-p_lead-2, 2]);
end

sgtitle('Problem 4(g): Lead Compensation Root Loci (z=1, varying p)');

% Open Control System Designer for each case
fprintf('\nOpening Control System Designer for each lead compensation case...\n');
for idx = 1:length(p_values)
    p_lead = p_values(idx);
    G_lead = tf([1 z_lead], conv([1 0 0], [1 p_lead]));
    fprintf('Opening designer for p = %d...\n', p_lead);
    controlSystemDesigner('rlocus', G_lead);
    pause(1); % Brief pause between windows
end

%% Part (h): Effect of moving pole
fprintf('\n--- Part (h): Effect of Moving Pole Location ---\n');
fprintf('As pole p moves RIGHT (decreases from 20 → 9 → 3):\n');
fprintf('- Pole approaches the zero at s = -1\n');
fprintf('- This is called "pole-zero cancellation" when p → z\n');
fprintf('- Root locus branches are "pulled" toward the pole\n');
fprintf('- Branches become more damped (move further left)\n');
fprintf('- Dominant poles have larger real parts → faster settling time\n');
fprintf('- When p is close to z, the compensator has less effect\n');
fprintf('\nEffect on response:\n');
fprintf('- p = 20: Pole far from zero, minimal cancellation, good lead effect\n');
fprintf('- p = 9: Moderate separation, moderate lead effect\n');
fprintf('- p = 3: Pole close to zero, near-cancellation, lead effect diminishes\n');
fprintf('- Optimal: Choose p > z with sufficient separation for desired response\n');

%% ========================================================================
%% PROBLEM 5: Motor Control (Previously completed)
%% ========================================================================
fprintf('\n========================================\n');
fprintf('PROBLEM 5: Motor Control System\n');
fprintf('========================================\n');

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
wn_min = 2.2 / (zeta_min * tr_spec);

fprintf('\nDesign Requirements:\n');
fprintf('Minimum damping ratio (zeta): %.4f\n', zeta_min);
fprintf('Minimum natural frequency (wn): %.2f rad/s\n', wn_min);

%% Part (a): Open sisotool with design requirements
fprintf('\n--- Part (a): SISO Design Tool ---\n');

controlSystemDesigner('rlocus', G5);

fprintf('Instructions for Part (a):\n');
fprintf('1. In Control System Designer window:\n');
fprintf('2. Right-click on root locus plot\n');
fprintf('3. Select "Design Requirements" -> "New..."\n');
fprintf('4. Add requirements:\n');
fprintf('   - Damping ratio >= %.4f (for Mp <= 5%%)\n', zeta_min);
fprintf('   - Natural frequency >= %.2f rad/s (for tr <= 0.5 ms)\n', wn_min);
fprintf('5. Shaded region shows acceptable pole locations\n');

%% Part (b): Select gain to meet specifications
fprintf('\n--- Part (b): Select Gain ---\n');
fprintf('Drag the pink square on root locus to select K\n');
fprintf('Place poles in shaded region\n');
fprintf('Note gain K and closed-loop poles at bottom of window\n');

% Example gain (replace with your selected value)
K_b = 50;
sys_cl_b = feedback(K_b * G5, 1);
poles_cl_b = pole(sys_cl_b);

fprintf('\nExample: K = %.2f\n', K_b);
fprintf('Closed-loop poles:\n');
disp(poles_cl_b);

%% Part (c): Step response and steady-state error
fprintf('\n--- Part (c): Step Response and Error ---\n');

figure('Name', 'Problem 5c - Step Response', 'NumberTitle', 'off');
step(sys_cl_b);
grid on;
title(sprintf('Closed-Loop Step Response (K = %.2f)', K_b));
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');

% Steady-state error
Kp_b = K_b * dcgain(G5);
ess_b = 1 / (1 + Kp_b);
ess_percent_b = ess_b * 100;

fprintf('Position constant Kp = %.4f\n', Kp_b);
fprintf('Steady-state error = %.4f (%.2f%%)\n', ess_b, ess_percent_b);

% Verify with FVT
T_b = feedback(K_b * G5, 1);
ess_fvt = 1 - dcgain(T_b);
fprintf('Steady-state error (FVT) = %.4f (%.2f%%)\n', ess_fvt, ess_fvt*100);

%% Part (d): Increase gain for 2% error
fprintf('\n--- Part (d): Gain for 2%% Error ---\n');

% Required gain
K_d = (1/0.02 - 1) / dcgain(G5);
fprintf('Required gain for 2%% error: K = %.2f\n', K_d);

sys_cl_d = feedback(K_d * G5, 1);
poles_cl_d = pole(sys_cl_d);

fprintf('Closed-loop poles:\n');
disp(poles_cl_d);

figure('Name', 'Problem 5d - Higher Gain Response', 'NumberTitle', 'off');
step(sys_cl_d);
grid on;
title(sprintf('Closed-Loop Step Response (K = %.2f, ess = 2%%)', K_d));
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');

% Check specifications
S_d = stepinfo(sys_cl_d);
fprintf('\nStep response characteristics:\n');
fprintf('  Rise time: %.4f ms\n', S_d.RiseTime * 1000);
fprintf('  Overshoot: %.2f%%\n', S_d.Overshoot);
fprintf('  Settling time: %.4f ms\n', S_d.SettlingTime * 1000);

if S_d.RiseTime > tr_spec
    fprintf('  WARNING: Rise time exceeds specification!\n');
end
if S_d.Overshoot > Mp_spec * 100
    fprintf('  WARNING: Overshoot exceeds 5%% specification!\n');
end

%% Part (e): Lag compensator exploration
fprintf('\n--- Part (e): Lag Compensator ---\n');

sisotool('rlocus', G5);

fprintf('Instructions:\n');
fprintf('1. In Control System Designer, go to Compensator tab\n');
fprintf('2. Right-click in Compensator Editor\n');
fprintf('3. Add "Real Pole" (close to origin)\n');
fprintf('4. Add "Real Zero" (slightly further from origin than pole)\n');
fprintf('5. For lag: place pole < zero (both near origin)\n');
fprintf('6. Observe how lag compensator affects branches:\n');
fprintf('   - Lag pulls branches SLIGHTLY LEFT\n');
fprintf('   - Main effect: increases DC gain (reduces steady-state error)\n');
fprintf('   - Minor effect on transient response if pole/zero close to origin\n');

%% Part (f): Minimum controller gain
fprintf('\n--- Part (f): Minimum Controller Gain ---\n');

% For C(s) = Kc*(s+z)/(s+p), DC gain = Kc*z/p
% For 2%% error: 1/(1 + (Kc*z/p)*dcgain(G)) = 0.02
min_gain_ratio = (1/0.02 - 1) / dcgain(G5);

fprintf('For 2%% steady-state error:\n');
fprintf('Required: Kc*z/p >= %.2f\n', min_gain_ratio);
fprintf('\nThis is the minimum DC gain of the compensator.\n');

%% Part (g): Design lag controller
fprintf('\n--- Part (g): Lag Controller Design ---\n');

% Example design (adjust based on sisotool results)
z_lag = 50;      % Zero location
p_lag = 5;       % Pole location (p < z for lag)
Kc = min_gain_ratio * p_lag / z_lag * 1.05;  % Add 5%% margin

fprintf('Example Lag Compensator:\n');
fprintf('  Zero (z) = %.2f\n', z_lag);
fprintf('  Pole (p) = %.2f\n', p_lag);
fprintf('  Kc = %.4f\n', Kc);
fprintf('  DC gain (Kc*z/p) = %.2f\n', Kc*z_lag/p_lag);

% Create compensator
C_lag = tf(Kc * [1, z_lag], [1, p_lag]);
G_comp = C_lag * G5;

fprintf('\nCompensated System:\n');
disp(G_comp);

% Closed-loop
sys_cl_lag = feedback(G_comp, 1);
poles_cl_lag = pole(sys_cl_lag);

fprintf('Closed-loop poles:\n');
disp(poles_cl_lag);

% Plot step response
figure('Name', 'Problem 5g - Lag Compensator Response', 'NumberTitle', 'off');
step(sys_cl_lag);
grid on;
title(sprintf('Step Response with Lag Compensator\n(Kc=%.4f, z=%.2f, p=%.2f)', Kc, z_lag, p_lag));
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');

% Performance metrics
S_lag = stepinfo(sys_cl_lag);
ess_lag = 1 - dcgain(sys_cl_lag);

fprintf('\nPerformance with lag compensator:\n');
fprintf('  Rise time: %.4f ms (spec: <= %.2f ms)\n', S_lag.RiseTime*1000, tr_spec*1000);
fprintf('  Overshoot: %.2f%% (spec: <= %.0f%%)\n', S_lag.Overshoot, Mp_spec*100);
fprintf('  Settling time: %.4f ms\n', S_lag.SettlingTime*1000);
fprintf('  Steady-state error: %.4f (%.2f%%)\n', ess_lag, ess_lag*100);

% Specification check
fprintf('\n** Specification Check **\n');
meets_tr = S_lag.RiseTime <= tr_spec;
meets_os = S_lag.Overshoot <= Mp_spec * 100;
meets_ess = ess_lag <= 0.02;

if meets_tr
    fprintf('✓ Rise time: PASS\n');
else
    fprintf('✗ Rise time: FAIL\n');
end
if meets_os
    fprintf('✓ Overshoot: PASS\n');
else
    fprintf('✗ Overshoot: FAIL\n');
end
if meets_ess
    fprintf('✓ Steady-state error: PASS\n');
else
    fprintf('✗ Steady-state error: FAIL\n');
end

% Root locus with compensator
figure('Name', 'Problem 5g - Root Locus with Lag', 'NumberTitle', 'off');
rlocus(G_comp);
grid on;
title('Root Locus with Lag Compensator');
xlabel('Real Axis');
ylabel('Imaginary Axis');
sgrid(zeta_min, wn_min);

%% Summary
fprintf('\n========================================\n');
fprintf('HOMEWORK COMPLETE\n');
fprintf('========================================\n');
fprintf('\nREMINDER:\n');
fprintf('1. Save all figures as image files\n');
fprintf('2. Include this code in your submission\n');
fprintf('3. Label all plots appropriately\n');
fprintf('4. For Problems 1-4, update the transfer functions in the code\n');
fprintf('5. Complete the MATLAB tutorial on Motor Speed Control\n');
fprintf('6. Adjust lag compensator parameters in sisotool for Problem 5g\n');
fprintf('\nGood luck!\n')