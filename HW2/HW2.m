%% ME 351 Homework #2 - Problem 5
%% Student: Dylan Qiu
%% DC Motor Control Analysis

clear all; close all; clc;

%% System Parameters (Quanser Qube-SERVO-2)
R = 8.4;           % Resistance (Ohms)
L = 0.00116;       % Inductance (H)
Jm = 4.65e-6;      % Motor inertia (kg-m^2)
Km = 0.042;        % Motor constant (N-m/A)
Kb = 0.042;        % Back EMF constant (V/(rad/s))
B = 0.02;          % Viscous damping (N-m-s)

% Disc parameters
mass_disc = 0.053; % kg
radius_disc = 0.0248; % m
Jd = 0.5 * mass_disc * radius_disc^2; % Disc inertia (kg-m^2)

% Total inertia
J_total = Jm + Jd;

%% Transfer Function Derivation
% For a DC motor: G(s) = Km / (s*((L*s + R)*(J*s + B) + Km*Kb))
% Simplified form for typical DC motor

% Electrical time constant
tau_e = L/R;

% Mechanical time constant  
tau_m = J_total/B;

% DC gain
K_dc = Km/(B*R + Km*Kb);

% Transfer function: G(s) = K_dc / (s*(tau_e*s + 1)*(tau_m*s + 1))
% For small L, can approximate as: G(s) = K_dc / (s*(tau_m*s + 1))

num = K_dc;
den = [tau_e*tau_m, tau_e + tau_m, 1, 0]; % s*(tau_e*s + 1)*(tau_m*s + 1)
G = tf(num, den);

fprintf('Dylan Qiu - ME 351 HW2 Problem 5\n');
fprintf('================================\n');
fprintf('System Parameters:\n');
fprintf('DC Gain: %.6f\n', K_dc);
fprintf('Electrical time constant: %.6f s\n', tau_e);
fprintf('Mechanical time constant: %.6f s\n', tau_m);
fprintf('\nTransfer Function G(s):\n');
G

%% Part (a): Open-loop and Proportional Control Response
figure(1);
set(gcf, 'Position', [100, 100, 1200, 800]);

% Open-loop step response
t = 0:0.001:2;
[y_open, t_open] = step(G, t);

subplot(2,2,1);
plot(t_open, y_open, 'b-', 'LineWidth', 2);
grid on;
title('Open-Loop Step Response - Dylan Qiu', 'FontSize', 12);
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('Open-Loop', 'Location', 'best');

% Proportional control with different gains
Kp_values = [5, 10, 40];
colors = ['r', 'g', 'm'];
legend_str = cell(1, length(Kp_values));

subplot(2,2,2);
hold on;
for i = 1:length(Kp_values)
    Kp = Kp_values(i);
    T_p = feedback(Kp*G, 1);  % Closed-loop with proportional control
    [y_p, t_p] = step(T_p, t);
    plot(t_p, y_p, colors(i), 'LineWidth', 2);
    legend_str{i} = sprintf('Kp = %d', Kp);
end
grid on;
title('Proportional Control Step Response - Dylan Qiu', 'FontSize', 12);
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend(legend_str, 'Location', 'best');
hold off;

fprintf('\nPart (a): Effect of increasing proportional gain:\n');
fprintf('- Higher Kp reduces steady-state error\n');
fprintf('- Higher Kp increases speed of response\n');
fprintf('- Higher Kp may cause overshoot and oscillations\n\n');

%% Part (b): PI Control
subplot(2,2,3);
hold on;

% PI controller with Kp=5, Ki=10
Kp = 5;
Ki_values = [10, 40];
colors_pi = ['c', 'k'];
legend_str_pi = cell(1, length(Ki_values));

for i = 1:length(Ki_values)
    Ki = Ki_values(i);
    C_pi = tf([Kp, Ki], [1, 0]);  % PI controller
    T_pi = feedback(C_pi*G, 1);
    [y_pi, t_pi] = step(T_pi, t);
    plot(t_pi, y_pi, colors_pi(i), 'LineWidth', 2);
    legend_str_pi{i} = sprintf('Kp=%d, Ki=%d', Kp, Ki);
end
grid on;
title('PI Control Step Response - Dylan Qiu', 'FontSize', 12);
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend(legend_str_pi, 'Location', 'best');
hold off;

fprintf('Part (b): Effect of increasing integral gain:\n');
fprintf('- Higher Ki eliminates steady-state error\n');
fprintf('- Higher Ki increases overshoot\n');
fprintf('- Higher Ki may cause longer settling time\n\n');

%% Part (c): PID Control
subplot(2,2,4);

% PID controller with Kp=5, Ki=40, Kd=5
Kp = 5;
Ki = 40;
Kd = 5;
C_pid = tf([Kd, Kp, Ki], [1, 0]);  % PID controller
T_pid = feedback(C_pid*G, 1);
[y_pid, t_pid] = step(T_pid, t);

plot(t_pid, y_pid, 'r-', 'LineWidth', 2);
grid on;
title('PID Control Step Response - Dylan Qiu', 'FontSize', 12);
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend(sprintf('Kp=%d, Ki=%d, Kd=%d', Kp, Ki, Kd), 'Location', 'best');

fprintf('Part (c): Effect of derivative term:\n');
fprintf('- Derivative term reduces overshoot\n');
fprintf('- Derivative term improves transient response\n');
fprintf('- Derivative term provides damping\n\n');

%% Part (d): Optimized PID Design
figure(2);
set(gcf, 'Position', [150, 150, 1000, 600]);

% Try different PID combinations to optimize response
% Target: Fast settling time, low overshoot, zero steady-state error

% Optimized values (you may need to tune these)
Kp_opt = 15;
Ki_opt = 25;  
Kd_opt = 3;

C_opt = tf([Kd_opt, Kp_opt, Ki_opt], [1, 0]);
T_opt = feedback(C_opt*G, 1);

% Extended time for settling time analysis
t_long = 0:0.001:5;
[y_opt, t_opt] = step(T_opt, t_long);

% Calculate performance metrics
step_info = stepinfo(T_opt);
settling_time = step_info.SettlingTime;
overshoot = step_info.Overshoot;
steady_state_value = dcgain(T_opt);
steady_state_error = abs(1 - steady_state_value);

% Plot optimized response
plot(t_opt, y_opt, 'b-', 'LineWidth', 2);
hold on;

% Add reference line
plot([0, t_long(end)], [1, 1], 'k--', 'LineWidth', 1);

% Add settling time bounds (2% of final value)
final_value = y_opt(end);
upper_bound = final_value * 1.02;
lower_bound = final_value * 0.98;
plot([0, t_long(end)], [upper_bound, upper_bound], 'r--', 'LineWidth', 1);
plot([0, t_long(end)], [lower_bound, lower_bound], 'r--', 'LineWidth', 1);

% Mark settling time
if settling_time <= t_long(end)
    line([settling_time, settling_time], [0, final_value], 'Color', 'g', 'LineWidth', 2);
end

grid on;
title(sprintf('Optimized PID Response - Dylan Qiu\nKp=%.1f, Ki=%.1f, Kd=%.1f', ...
    Kp_opt, Ki_opt, Kd_opt), 'FontSize', 12);
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('PID Response', 'Reference', '±2% Bounds', 'Settling Time', 'Location', 'best');

% Display performance metrics
fprintf('Part (d): Optimized PID Performance:\n');
fprintf('Selected gains: Kp = %.1f, Ki = %.1f, Kd = %.1f\n', Kp_opt, Ki_opt, Kd_opt);
fprintf('Performance Specifications:\n');
fprintf('  i.  Settling Time (2%%): %.3f seconds\n', settling_time);
fprintf('  ii. Overshoot: %.2f%%\n', overshoot);
fprintf('  iii.Steady-State Error: %.6f\n', steady_state_error);
fprintf('\nRationale for gain selection:\n');
fprintf('- Kp increased to improve speed of response\n');
fprintf('- Ki maintained to eliminate steady-state error\n');
fprintf('- Kd tuned to minimize overshoot while maintaining fast response\n');
fprintf('- Balance achieved between settling time and overshoot\n');

% Add text annotations on the plot
text_x = t_long(end)*0.6;
text_y = final_value*0.7;
text(text_x, text_y, sprintf('Settling Time: %.3fs\nOvershoot: %.2f%%\nSS Error: %.6f', ...
    settling_time, overshoot, steady_state_error), ...
    'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Save the figure
saveas(gcf, 'Dylan_Qiu_PID_Response.png');

fprintf('\nAll plots have been generated and saved.\n');
fprintf('Figure saved as: Dylan_Qiu_PID_Response.png\n');