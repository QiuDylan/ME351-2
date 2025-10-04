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
K = 0.01;
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


fprintf('PROBLEM 4d\n');
% Transfer function: G1(s) = s^2 + K(s+1)
K = 0.01;
J = 1 ; %kg-m^2
s = tf('s');
G4d = K * (s + 1) / s^2 ;
%controlSystemDesigner('rlocus', G4d)

fprintf('\nPart (g): Lead Compensation with z=1, varying p \n');

z_lead = 1;
p_values = [20, 9, 3];

figure('Name', 'Problem 4g - Lead Compensation Comparison\n', 'NumberTitle', 'off');

for idx = 1:length(p_values)
    p_lead = p_values(idx);
    
    fprintf('\n(g.%d) p = %d:\n', idx, p_lead);
    
    % Transfer function
    sg = tf('s');
    G_lead = (sg + 1) / (sg^3 + p_lead * sg^2);
    G_lead
    
    % Plot
    subplot(1, 3, idx);
    %rlocus(G_lead);
    %grid on;
    title(sprintf('p = %d', p_lead));
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    %sgrid;
    %axis equal;
    %xlim([-p_lead-2, 2]);
end

