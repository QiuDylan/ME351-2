% ME351 In-class Motor Problem Matlab Code
% Root locus and Lead/lag example
% Dylan Qiu, ME '27

clear all; close all; clc;

K = 1;
tau = 1;

% Input transfer
figure(1);
num = 1;
den = [tau 1 0];
sys = tf(num,den);

step(sys)
title('Open-loop response');
%cltf = feedback(sys)
rlocus(sys)

% In factored form 
