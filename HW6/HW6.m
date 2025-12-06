% Dylan Qiu ME 351 Homework #6 
% Clear workspace
clear all; close all; clc;

%% Problem 1: G(s) = 1/[s(s^2 + 0.8s + 1)]
fprintf('=== PROBLEM 1 ===\n');

% Define transfer function
num1 = 1;
den1 = conv([1 0], [1 0.8 1]);  % s * (s^2 + 0.8s + 1)
G1 = tf(num1, den1);

% Plot Nyquist diagram
figure(1);
nyquist1(G1);
grid on;
title('Problem 1: Nyquist Diagram for G(s) = 1/[s(s^2 + 0.8s + 1)]');

% Stability analysis
fprintf('\nProblem 1 Analysis:\n');
fprintf('Open-loop poles: ');
poles1 = pole(G1);
disp(poles1);
fprintf('Number of open-loop poles in RHP (P): %d\n', sum(real(poles1) > 0));
fprintf('The system has 1 pole at origin and 2 complex conjugate poles in LHP\n');
fprintf('For stability: N (encirclements of -1) = -P = 0\n');
fprintf('The Nyquist plot does NOT encircle -1, so N = 0\n');
fprintf('Therefore: Z = N + P = 0 + 0 = 0 (STABLE closed-loop system)\n\n');

%% Problem 2: G(s) = (s^2 + 2s + 1)/(s^3 + 0.2s^2 + s + 1)
fprintf('=== PROBLEM 2 ===\n');

% Define transfer function
num2 = [1 2 1];
den2 = [1 0.2 1 1];
G2 = tf(num2, den2);

% Plot Nyquist diagram
figure(2);
nyquist1(G2);
grid on;
title('Problem 2: Nyquist Diagram for G(s) = (s^2+2s+1)/(s^3+0.2s^2+s+1)');

% Stability analysis
fprintf('\nProblem 2 Analysis:\n');
fprintf('Open-loop poles: ');
poles2 = pole(G2);
disp(poles2);
P2 = sum(real(poles2) > 0);
fprintf('Number of open-loop poles in RHP (P): %d\n', P2);
fprintf('Check if Nyquist plot encircles -1 point\n');
fprintf('From the plot, count encirclements (clockwise positive)\n');
fprintf('For stability: Z = N + P should equal 0\n\n');

%% Problem 3: Inverted Pendulum
fprintf('=== PROBLEM 3: INVERTED PENDULUM ===\n');

% Given parameters
Jp = 1.20e-4;      % kg-m^2
Mp = 0.0270;       % kg
r = 0.0826;        % m
lp = 0.153;        % m
Jeq = 1.23e-4;     % kg-m^2
R = 3.30;          % Ohms
Km = 0.028;        % Motor constant
g = 9.81;          % m/s^2

% Define transfer function
num = [Km*r*Mp*lp 0 0];
den = conv([1 0], conv([R*Jeq Km^2], [Jp 0 -g*Mp*lp]));
invpend = tf(num, den);

fprintf('\nInverted Pendulum Transfer Function:\n');
invpend

% Part (a): Uncompensated Bode and Nyquist plots
fprintf('\n--- Part (a): Uncompensated Plots ---\n');

figure(3);
subplot(2,1,1);
bode(invpend);
grid on;
title('Part (a): Uncompensated Bode Plot');

subplot(2,1,2);
nyquist1(invpend);
grid on;
title('Part (a): Uncompensated Nyquist Plot');

% Part (b): Comment on encirclements
fprintf('\n--- Part (b): Stability Analysis ---\n');
poles_invpend = pole(invpend);
fprintf('Open-loop poles:\n');
disp(poles_invpend);
P_invpend = sum(real(poles_invpend) > 0);
fprintf('Number of open-loop RHP poles (P): %d\n', P_invpend);

%% Part (c): Add pole at origin
fprintf('\n--- Part (c): Add Pole at Origin ---\n');

C_pole = tf(1, [1 0]);  % 1 / s
GC_c = invpend * C_pole;

figure(4);
nyquist1(GC_c);
grid on;
title('Part (c): Nyquist with Additional Pole at Origin');


%% Part (d): Lead-Lag Controller Design
fprintf('\n--- Part (d): Lead-Lag Controller Design ---\n');

% Lag controller parameters
z_lag = 3;
p_lag = 0;  

% Lead controller parameters
z_lead = 4;
w_max = 20;  % rad/s 

% Calculate p_lead for max phase at 20 rad/s
% For lead compensator: w_max = sqrt(z_lead * p_lead)
p_lead = (w_max^2) / z_lead;

fprintf('Lead controller parameters:\n');
fprintf('  z_lead = %.2f rad/s\n', z_lead);
fprintf('  p_lead = %.2f rad/s\n', p_lead);
fprintf('  Maximum phase occurs at w = %.2f rad/s\n', w_max);

% Maximum phase lead
phi_max = asin((p_lead - z_lead)/(p_lead + z_lead)) * 180/pi;
fprintf('\n Maximum phase lead Phi = %.2f degrees\n', phi_max);

% Define controller 
C_lag = tf([1 z_lag], [1 p_lag]);
C_lead = tf([1 z_lead], [1 p_lead]);
C_nok = C_lag * C_lead;

% Open-loop with controller (K=1)
GC_d = invpend * C_nok;

figure(5);
subplot(2,1,1);
bode(GC_d);
grid on;
title('Part (d): Bode Plot with Lead-Lag (K=1)');

subplot(2,1,2);
nyquist1(GC_d);
grid on;
title('Part (d): Nyquist with Lead-Lag (K=1)');

%% Part (e): Adjust gain for PM > 50 degrees
fprintf('\n--- Part (e): Gain Selection for PM > 50 degrees ---\n');

% Find gain for desired phase margin
controlSystemDesigner(GC_d);
%% Part (f): Plot with selected gain
fprintf('\n--- Part (f): Plots with Selected Gain ---\n');
best_K = 34.7; %from control system designer
C_final = best_K * C_nok;
GC_final = best_K * GC_d;

[Gm_final, Pm_final, Wcg_final, Wcp_final] = margin(GC_final);

fprintf('Final System Margins:\n');
fprintf('  Gain Margin: %.2f dB at %.2f rad/s\n', 20*log10(Gm_final), Wcg_final);
fprintf('  Phase Margin: %.2f degrees at %.2f rad/s\n\n', Pm_final, Wcp_final);

figure(6);
subplot(2,1,1);
margin(GC_final);
grid on;
title(sprintf('Part (f): Bode Plot with K=%.4f (PM=%.2f deg)', best_K, Pm_final));

subplot(2,1,2);
nyquist1(GC_final);
grid on;
title(sprintf('Part (f): Nyquist with K=%.4f', best_K));

%% Part (g): Closed-loop disturbance impulse response
fprintf('\n--- Part (g): Closed-Loop Disturbance Response ---\n');


% Y(s)/D(s) = G(s)/(1 + G(s)C(s))
T_disturbance = feedback(invpend, C_final);

figure(7);
impulse(T_disturbance);
grid on;
title('Part (g): Closed-Loop Impulse Response to Disturbance');
xlabel('Time (seconds)');
ylabel('Pendulum Angle (radians)');

