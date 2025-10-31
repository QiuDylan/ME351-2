% ME 351 - Homework 4: Bode Plots
% Problems 3 and 4.

% Problem 3
figure(1);
% 3a: G(s) = 10/s

s = tf('s');
G3a = 10/s;
subplot(2, 2, 1);
bode(G3a, bodeoptions);
title('Problem 3a');
grid on;

% 3b: G(s) = s^2
G3b = s^2;
subplot(2, 2, 2);
bode(G3b, bodeoptions);
title('Problem 3b');
grid on;

% 3c: G(s) = 5 / (s + 5)
G3c = 5 / (s + 5);
subplot(2, 2, 3);
bode(G3c, bodeoptions);
title('Problem 3c');
grid on;

% 3d: G(s) = 10s / (s + 10)
G3d = (10*s) / (s + 10);
subplot(2, 2, 4);
bode(G3d, bodeoptions);
title('Problem 3d');
grid on;

% Problem 4
figure(2);
% 4a: G(s) = 100 / (s^2 + 2s + 100)
G4a = 100 / (s^2 + 2*s + 100);
subplot(3, 1, 1);
bode(G4a, bodeoptions);
title('Problem 4a');
grid on;

% 4b: G(s) = (s + 5) / (s^2 + 2s + 100)
G4b = (s + 5) / (s^2 + 2*s + 100);
subplot(3, 1, 2);
bode(G4b, bodeoptions);
title('Problem 4b');
grid on;

% 4c: G(s) = (s^2 + 2s + 100) / ((s + 1) * (5s + 500))
G4c = (s^2 + 2*s + 100) / ((s + 1) * (5*s + 500));
subplot(3, 1, 3);
bode(G4c, bodeoptions);
title('Problem 4c');
grid on;
