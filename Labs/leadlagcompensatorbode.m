% Dylan Qiu 
clc; clear; close all;

% lead compensator

figure(1);
s = tf('s');


G = (200)/ ((s+0.1)*(s^2+5*s+100));
feedback(G,1);
step(G)
figure(2);
%{
z_lag = 10;
p_lag = 0.1;
z_lead = 5;
p_lead = 10;
K1 = 10;
K2 = 5/4;
C_lag = K1 * (s/z_lag+1)/(s/p_lag+1); %lag compensator 
C_lead = K2 * (s+z_lead)/(s+p_lead); %lead compensator 
G_comp = G * C_lead;
%}
bode(G, 'y',bodeoptions); 
legend on; 