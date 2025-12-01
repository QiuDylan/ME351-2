%% Dylan Qiu, ME '27, ME351, Prof. Baglione
%% Midterm 1: Take Home portion

% Question 1
% code borrowed from https://ctms.engin.umich.edu/CTMS/index.php?aux=Extras_Ess

Kp = (7 / 1.2)^2;
z = 1;
p = z * 0.01 * Kp / 7;

s = tf('s');
G = (Kp)/(s*(s+7)); % open loop tf
C = (s+z)/(s+p); %lag compensator 
%C = 1; % run this line to get only Kp
T = feedback(G * C,1);
t = 0:0.1:4;
u = t;
[y,t,x] = lsim(T,u,t);
figure(1);
plot(t,y,'y',t,u,'m')
xlabel('Time (sec)')
ylabel('Amplitude')
title('Kp and lag compensator, Input-purple, Output-yellow')

% Question 2

K = 40; % or 36 in part b
a = 8;
s = tf('s');
G2 = 1 /(s^2 + (a - 4) * s - 4 * a + K);
sys2 = K * G2;
figure(2);
step(sys2)
stats = stepinfo(sys2);
fprintf('Tr = %4.2f\n', stats.RiseTime);
fprintf('Overshoot = %4.2f\n', stats.Overshoot);

% Question 3

z_values = [2, 1, 4, 20];

%fprintf('Problem 3a\n');

for idx = 1:length(z_values)
    z_lead = z_values(idx);
    
    fprintf('\n(g.%d) p = %d:\n', idx, z_lead);
    
    % Transfer function
    sg = tf('s');
    G3 = (sg + z_lead) / (sg^2 + 4 * sg + 8);
    G3;
    figure(3);
    subplot(1, length(z_values), idx);
    rlocus(G3);
    title(sprintf('zeroes = %d', z_lead));
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    axis equal;
end
%fprintf('Problem 3b\n');
results = zeros(length(z_values), 4);

for i = 1:length(z_values)
    z = z_values(i);
   
    a = 1;
    b = 8 - 4 * z;
    c = -16;
    K = (-b + sqrt(b^2 - 4 * a * c)) / (2 * a); % Solve for K
    
    den = s^2 + 4 * s + 8;
    P = (s + z) / den;
    CLTF = feedback(K * P, 1);
   
    S = stepinfo(CLTF);
   
    ess = 8 / (8 + K * z);
    yss = K * z / (8 + K * z);
    results(i, :) = [z, K, S.SettlingTime, S.Overshoot];
    fprintf('| z = %-4.0f | K = %-5.2f | Ts = %-5.3f s | y_ss = %-8.3f | Overshoot = %-8.3f%% | e_ss = %-4.3f |\n', ...
            z, K, S.SettlingTime, yss, S.Overshoot, ess);
end
