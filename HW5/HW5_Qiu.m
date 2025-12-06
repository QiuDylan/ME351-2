% Dylan Qiu ME 351 HW 5
clc; clear; close all;

% part a
%subplot(2,2,1);
figure(1);
s = tf('s');
num = 26;
den = [1/10, 1, 0];
G = tf(num,den);
margin(G);
grid on;
bode(G);
title('Bode Plot - Dylan Qiu');

%part d 
figure(2);
%subplot(2,2,2);
kp = 100 / 26; % new compensator gain
G_comp = tf(kp*num,den);
bode(G_comp, G);
legend on;
title('Bode Plot with Compensator - Dylan Qiu');

% Part e
figure(3);
G_step = feedback(kp * G,1);
step(G_step)

% Part f
figure(4);
%C = (s+z)/(s+p); %lag compensator 
C = 1; % run this line to get only Kp
T = feedback(G * kp,1);
t = 0:0.1:4;
u = t;
[y,t,x] = lsim(T,u,t);
plot(t,y,'y',t,u,'m')
xlabel('Time (sec)')
ylabel('Amplitude')
grid on;
title('Ramp Response, Input-purple, Output-yellow - Dylan Qiu')
%title('Ramp Response with Kp - Dylan Qiu');

% Part G
figure(5);
OS_max = 20;  % percent
zeta_min = sqrt((log(OS_max/100))^2 / (pi^2 + (log(OS_max/100))^2));
pm = 100 * zeta_min; 

[~, ~, ~, Wcp_uncomp] = margin(G);
Wcp_comp = 2 * Wcp_uncomp; % new crossover frequency
[mag_wcp, phase_wcp] = bode(G, Wcp_comp);
fprintf('Desired compensated crossover freq: %.4f rad/s\n', Wcp_comp);

p_max = pm - (180 + phase_wcp) + 5;
phi_max = p_max * pi/180;
alpha = (1 - sin(phi_max))/(1 + sin(phi_max));
z = Wcp_comp * sqrt(alpha);
p = z / alpha;

fprintf('mag_wcp %.4f, phase_wcp %.4f\n', mag_wcp, phase_wcp);
Kc = 1 / (mag_wcp * sqrt(alpha));
fprintf('Compensator gain Kc: %.4f\n', Kc);
Comp = tf(Kc*[1,z], [1,p]);
fprintf('\nLead Compensator: Gc(s) = %.4f * (s + %.4f)/(s + %.4f)\n', Kc, z, p);
G_lead = G * Comp ;
bode(G_lead,G,G_comp);
legend on;
title('Bode Plot Comparison - Dylan Qiu');
[Gm_comp, Pm_comp, Wcg_comp, Wcp_comp_actual] = margin(G_lead);

fprintf('Phase Margin: %.2f degrees\n', Pm_comp);
fprintf('Gain Margin: %.2f dB\n', 20*log10(Gm_comp));
fprintf('Crossover Frequency: %.4f rad/s\n', Wcp_comp_actual);
G_step_comp= feedback(G* Comp,1);
figure(6);
step(G_step, 'r', G_step_comp, 'g');
legend on;
grid on;
title('Lead Compensated Step Response - Dylan Qiu');

figure(7);
%C = (s+z)/(s+p); %lag compensator 
%C = 1; % run this line to get only Kp
%T = feedback(G * kp,1);
t1 = 0:0.1:4;
u1 = t1;
T = feedback(G_lead,1);
[y1,t1,x1] = lsim(T,u1,t1);
plot(t1,y1,'y',t1,u1,'m')
xlabel('Time (sec)')
ylabel('Amplitude')
grid on;
title('Ramp Response, Input-purple, Output-yellow - Dylan Qiu')
%title('Ramp Response with Kp - Dylan Qiu');

ess_comp = 1 / (Kc * kp * 26);
fprintf('Steady-state error to ramp: %.6f (Required: < 0.01)\n', ess_comp);
fprintf('SS error spec met: %s\n', string(ess_comp < 0.01));

fprintf('Bandwidth: %.4f rad/s\n', bandwidth(T));