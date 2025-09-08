% ME 351 - DC Motor Step Response Simulation
% Quanser Qube-SERVO-2 Open-Loop Response
clear all; close all; clc;

%% System Parameters 
R = 8.4;        % Armature resistance (Ohms)
L = 0.00116;    % Armature inductance (H)
J_m = 4.65e-6;  % Motor inertia (kg-m^2)
K_m = 0.042;    % Motor torque constant (N-m/A)
K_b = 0.042;    % Back-emf constant (V/(rad/s))

% Disc properties
m_disc = 0.053;   % Disc mass (kg)
r_disc = 0.0248;  % Disc radius (m)
J_d = 0.5 * m_disc * r_disc^2;  % Disc inertia
J = J_m + J_d;    % Total inertia

fprintf('Disc inertia J_d = %.6e kg-m^2\n', J_d);
fprintf('Total inertia J = %.6e kg-m^2\n', J);

%% Three different damping values
B_values = [0, 0.0005, 0.005];  % N-m-s
colors = {'b', 'r', 'g'};
line_styles = {'-', '-', '-'};

%% Create figure for plotting
figure('Position', [100, 100, 800, 600]);
hold on;
grid on;

%% Simulate for each damping value
for i = 1:length(B_values)
    B = B_values(i);
    
    % Transfer function: Omega(s)/V(s) = K_m / ((Ls + R)(Js + B) + K_m*K_b)
    
    num = K_m;  % Numerator
    den = [L*J, (L*B + R*J), (R*B + K_m*K_b)];  % Denominator coefficients
    
    % Create transfer function
    sys = tf(num, den);
    
    % Generate step response (1V input)
    t = 0:0.001:0.5;  % Time vector from 0 to 0.5 seconds
    [y, t_out] = step(sys, t);
    
    % Calculate steady-state value analytically
    omega_ss = K_m / (R*B + K_m*K_b);
    
    % Plot the response
    plot(t_out, y, 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', sprintf('B = %.4f, \\omega_{ss} = %.2f rad/s', B, omega_ss));
    
    % Steady-state value
    plot([0, 0.5], [omega_ss, omega_ss], '--', 'Color', colors{i}, 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Print results 
    fprintf('\n----- B = %.4f N-m-s -----\n', B);
    fprintf('Steady-state angular velocity (analytical): %.2f rad/s\n', omega_ss);
    fprintf('Steady-state angular velocity (simulated): %.2f rad/s\n', y(end));
    
    % Display transfer function
    fprintf('Transfer function:\n');
    disp(sys);
end

%% Format the plot
xlabel('Time (s)', 'FontSize', 12);
ylabel('Angular Velocity \omega (rad/s)', 'FontSize', 12);
title('DC Motor Open-Loop Unit Step Response - Quanser Qube-SERVO-2', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
xlim([0, 0.5]);
ylim([0, 30]);

%{ 
Add text box with system parameters
str = sprintf(['System Parameters:\n' ...
               'R = %.1f \\Omega\n' ...
               'L = %.4f H\n' ...
               'J = %.6e kg-m^2\n' ...
               'K_m = %.3f N-m/A\n' ...
               'K_b = %.3f V/(rad/s)'], ...
               R, L, J, K_m, K_b);
annotation('textbox', [0.65, 0.2, 0.25, 0.25], ...
           'String', str, ...
           'FontSize', 9, ...
           'BackgroundColor', 'white', ...
           'EdgeColor', 'black');
%}

%% Part (e): Compare steady-state values
fprintf('Part (e): Comparison of steady-state angular velocities\n');
for i = 1:length(B_values)
    B = B_values(i);
    omega_ss_analytical = K_m / (R*B + K_m*K_b);
    fprintf('B = %.4f N-m-s: omega_ss = %.2f rad/s\n', B, omega_ss_analytical);
end

% saveas(gcf, 'DC_Motor_Step_Response.png');
% saveas(gcf, 'DC_Motor_Step_Response.fig');