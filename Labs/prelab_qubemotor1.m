%% ME352 Pre-Lab Assignment Parts 5 & 6 - MATLAB Calculations
% DC Motor Parameter Analysis
% Authors: [Your Names]
% Date: [Current Date]

clear all; clc; close all;

%% Given Parameters from Quanser Manual
fprintf('=== DC Motor Parameters ===\n');
R = 8.4;                    % Resistance [Ohm]
L = 0.00116;               % Inductance [H]
J_m = 4.65e-6;             % Motor inertia [kg-m^2]
K_m = 0.042;               % Motor torque constant [N-m/A]
K_b = 0.042;               % Back EMF constant [V/(rad/s)]
B = 0;                     % Viscous damping (negligible)

% Disc parameters
disc_mass = 0.053;         % Disc mass [kg]
disc_radius = 0.0248;      % Disc radius [m]

fprintf('R = %.1f Ohm\n', R);
fprintf('L = %.5f H\n', L);
fprintf('J_m = %.2e kg-m^2\n', J_m);
fprintf('K_m = %.3f N-m/A\n', K_m);
fprintf('K_b = %.3f V/(rad/s)\n', K_b);
fprintf('Disc mass = %.3f kg\n', disc_mass);
fprintf('Disc radius = %.4f m\n', disc_radius);

%% Part 5: Calculate disc inertia and total inertia
fprintf('\n=== PART 5: Inertia and Pole Calculations ===\n');

% Calculate disc inertia (solid disc: J = 1/2 * m * r^2)
J_d = 0.5 * disc_mass * disc_radius^2;
fprintf('Disc inertia J_d = %.2e kg-m^2\n', J_d);

% Calculate total inertia
J = J_m + J_d;
fprintf('Total inertia J = J_m + J_d = %.2e kg-m^2\n', J);

% Second-order characteristic equation: JL*s^2 + JR*s + K_m*K_b = 0
% Coefficients of the characteristic polynomial
a2 = J * L;                % s^2 coefficient
a1 = J * R;                % s^1 coefficient  
a0 = K_m * K_b;            % s^0 coefficient

fprintf('\nCharacteristic equation: %.2es^2 + %.2es + %.2e = 0\n', a2, a1, a0);

% Calculate poles using quadratic formula
discriminant = a1^2 - 4*a2*a0;
if discriminant >= 0
    pole1 = (-a1 + sqrt(discriminant)) / (2*a2);
    pole2 = (-a1 - sqrt(discriminant)) / (2*a2);
    fprintf('Poles are real:\n');
    fprintf('  pole1 = %.1f rad/s\n', pole1);
    fprintf('  pole2 = %.1f rad/s\n', pole2);
else
    real_part = -a1 / (2*a2);
    imag_part = sqrt(-discriminant) / (2*a2);
    fprintf('Poles are complex:\n');
    fprintf('  pole1 = %.1f + j%.1f rad/s\n', real_part, imag_part);
    fprintf('  pole2 = %.1f - j%.1f rad/s\n', real_part, imag_part);
end

% Alternative calculation using MATLAB roots function
coeffs = [a2, a1, a0];
poles_matlab = roots(coeffs);
fprintf('\nVerification using MATLAB roots():\n');
for i = 1:length(poles_matlab)
    if imag(poles_matlab(i)) == 0
        fprintf('  pole%d = %.1f rad/s\n', i, real(poles_matlab(i)));
    else
        fprintf('  pole%d = %.1f %+.1fj rad/s\n', i, real(poles_matlab(i)), imag(poles_matlab(i)));
    end
end

% Analyze pole dominance
fprintf('\nPole Analysis:\n');
if discriminant < 0
    fprintf('System is underdamped (complex poles)\n');
    fprintf('Natural frequency ωn = %.1f rad/s\n', sqrt(a0/a2));
    fprintf('Damping ratio ζ = %.3f\n', a1/(2*sqrt(a2*a0)));
    
    % Compare electrical vs mechanical time constants
    tau_electrical = L/R;
    tau_mechanical = J*R/(K_m*K_b);
    fprintf('Electrical time constant τ_e = L/R = %.2e s\n', tau_electrical);
    fprintf('Mechanical time constant τ_m = JR/(KmKb) = %.2e s\n', tau_mechanical);
    
    if tau_electrical < tau_mechanical
        fprintf('Electrical dynamics are faster → electrical pole dominates\n');
    else
        fprintf('Mechanical dynamics are faster → mechanical pole dominates\n');
    end
end

%% Part 6: First-order analysis (L = 0)
fprintf('\n=== PART 6: First-Order Analysis (L = 0) ===\n');

% Time constant and DC gain for first-order system
tau = (J * R) / (K_m * K_b);
K_dc = 1 / K_b;

fprintf('Time constant τ = JR/(K_m*K_b) = %.4f s\n', tau);
fprintf('DC gain K = 1/K_b = %.2f (rad/s)/V\n', K_dc);

% First-order pole location
pole_first_order = -1/tau;
fprintf('First-order pole = %.1f rad/s\n', pole_first_order);

% Physical interpretation
fprintf('\nPhysical Interpretation:\n');
fprintf('• Time to reach 63%% of final value: τ = %.3f s\n', tau);
fprintf('• Time to reach 95%% of final value: 3τ = %.3f s\n', 3*tau);
fprintf('• Steady-state speed per volt: %.1f rad/s/V\n', K_dc);
fprintf('• System bandwidth (approx): 1/τ = %.1f rad/s\n', 1/tau);

%% Create transfer functions and compare responses
fprintf('\n=== Transfer Function Comparison ===\n');

% Second-order transfer function
num_2nd = K_m;
den_2nd = [a2, a1, a0];
G_2nd = tf(num_2nd, den_2nd);

% First-order transfer function  
num_1st = K_dc;
den_1st = [tau, 1];
G_1st = tf(num_1st, den_1st);

fprintf('Second-order TF:\n');
G_2nd

fprintf('First-order TF:\n');
G_1st

%% Plot step responses for comparison
figure(1);
t = 0:0.001:0.5;  % Time vector for 0.5 seconds

% Step responses
[y_2nd, t_2nd] = step(G_2nd, t);
[y_1st, t_1st] = step(G_1st, t);

plot(t_2nd, y_2nd, 'b-', 'LineWidth', 2, 'DisplayName', '2nd Order (with L)');
hold on;
plot(t_1st, y_1st, 'r--', 'LineWidth', 2, 'DisplayName', '1st Order (L=0)');
grid on;
xlabel('Time [s]');
ylabel('Angular Velocity [rad/s]');
title('Step Response Comparison: 2nd Order vs 1st Order Models');
legend('Location', 'southeast');

% Add time constant markers
line([tau tau], [0 K_dc*0.632], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
text(tau+0.01, K_dc*0.632, sprintf('τ = %.3f s\n(63%% of final)', tau), 'FontSize', 10);

%% Summary table
fprintf('\n=== SUMMARY TABLE ===\n');
fprintf('Parameter                    | Value              | Units\n');
fprintf('-----------------------------|--------------------|-----------\n');
fprintf('Total inertia J              | %.2e        | kg-m^2\n', J);
fprintf('Time constant τ              | %.4f           | s\n', tau);
fprintf('DC gain K                    | %.2f            | (rad/s)/V\n', K_dc);
fprintf('First-order pole             | %.1f            | rad/s\n', pole_first_order);
if discriminant < 0
    fprintf('Second-order poles (complex) | %.1f±j%.1f      | rad/s\n', real(poles_matlab(1)), abs(imag(poles_matlab(1))));
end

fprintf('\nConclusion: The first-order approximation (L=0) provides a simpler\n');
fprintf('model that captures the essential motor dynamics for control purposes.\n');