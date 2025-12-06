% Dylan Qiu In class freq response workshop

clc; clear; close all;

% Part 1 
figure(1);
s = tf('s');

num = 1;
den = [1, 1, 0];
G = tf(num,den);
z_lead = 2;
p_lead = 10;
K2 = 66;
%C_lag = K1 * (s/z_lag+1)/(s/p_lag+1); %lag compensator 
C_lead = K2 * (s+z_lead)/(s+p_lead); %lead compensator 
G_comp = G * C_lead;

bode(G,'b', G_comp, 'r'); 
legend on; 

G_step= feedback(G_comp,1);
figure(2);
step(G_step, 'r');
legend on;
grid on;

% Part 2 
%{
figure(3);
t1 = 0:0.1:4;
u1 = t1;
T = G_step;
[y1,t1,x1] = lsim(T,u1,t1);
plot(t1,y1,'y',t1,u1,'m')
xlabel('Time (sec)')
ylabel('Amplitude')
grid on;
title('Ramp Response, Input-purple, Output-yellow - Dylan Qiu')
%}
K1 = 66*8;
z_lag = 0.0008;
p_lag = 0.0001;
C_lag = K1 *(z_lag/p_lag) * (s/z_lag+1)/(s/p_lag+1); %lag compensator 

figure(3);
t1 = 0:0.1:4;
u1 = t1;
T = feedback(G * C_lag,1);
[y1,t1,x1] = lsim(T,u1,t1);
plot(t1,y1,'y',t1,u1,'m')
xlabel('Time (sec)')
ylabel('Amplitude')
grid on;
title('Ramp Response, Input-purple, Output-yellow - Dylan Qiu')
ess = u1(end)-y1(end);
fprintf('final value, %.4f', ess);

K6 = 1/66;
G6_lag = G * C_lag * K6;
G6_lead = G * C_lead * K6;
G1_step = feedback(G,1);
G6lag_step= feedback(G6_lag,1);
G6lead_step= feedback(G6_lead,1);
%G6lead_lag_step = feedback(G6_lead*G6_lag,1);
figure(4);
step(G1_step,'y', G6lag_step, 'r', G6lead_step,'b');
legend on;
grid on;

% part 3

figure(5);

s3 = tf('s');
num3 = 10;
den3 = [1, 5, 4];
G3 = tf(num3,den3);
z3_lead = 1;
p3_lead = 1000;
K3 = 40;
C3_lead = K3 * (s/z3_lead+1)/(s/p3_lead+1); %lead compensator 
G3_comp = G3 * C3_lead;
bode(G3_comp, G3);

figure(6);

G3_step = feedback(G3_comp,1);
step(G3_step)
